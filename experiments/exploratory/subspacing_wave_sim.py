"""Truth-side geometry, numerical quadrature and tiled strict-order oracle.

No truth-side object belongs in the observer. Quadrature is convergent numerical
integration, not symbolic exactness. The ellipse area is 2*pi/sqrt(det K);
integrating its v width gives pi/sqrt(det K).
"""

from __future__ import annotations

import os

import numpy as np
from numpy.polynomial import legendre, polynomial

# Sandboxed Windows cache writes can stall tempfile; CUDA binaries may stay in RAM.
os.environ.setdefault("CUPY_CACHE_IN_MEMORY", "1")


class PanelPrimitive:
    """Composite Gauss polynomial antiderivative, including the local panel."""

    def __init__(self, function, panels, degree=8):
        self.panels = int(panels)
        nodes, weights = legendre.leggauss(degree)
        self.width = 1 / self.panels
        starts = np.arange(self.panels) * self.width - 0.5
        values = function(starts[:, None] + (nodes + 1) * self.width / 2)
        coeff = legendre.legfit(nodes, values.T, degree - 1)
        integ = legendre.legint(coeff, axis=0) * self.width / 2
        # One small basis conversion replaces a Python loop over every panel.
        conversion = np.zeros((degree + 1, degree + 1))
        for j in range(degree + 1):
            column = legendre.leg2poly(np.eye(degree + 1)[j])
            conversion[: len(column), j] = column
        self.coeff = integ.T @ conversion.T
        self.coeff[:, 0] -= polynomial.polyval(-1, self.coeff.T)
        totals = values @ weights * self.width / 2
        self.prefix = np.r_[0, np.cumsum(totals)]

    def __call__(self, u):
        u = np.asarray(u)
        if np.any((u < -0.5) | (u > 0.5)):
            raise ValueError("primitive domain is [-1/2,1/2]")
        flat = u.ravel()
        result = np.empty(flat.size, dtype=self.coeff.dtype)
        for start in range(0, flat.size, 262144):
            at = flat[start : start + 262144]
            scaled = (at + 0.5) * self.panels
            idx = np.minimum(scaled.astype(np.int64), self.panels - 1)
            x = 2 * (scaled - idx) - 1
            val = self.coeff[idx, -1]
            for j in range(self.coeff.shape[1] - 2, -1, -1):
                val = val * x + self.coeff[idx, j]
            result[start : start + len(at)] = self.prefix[idx] + val
        return result.reshape(u.shape)

    def gpu(self, u):
        """Fused Horner evaluation; caller synchronizes for honest timing."""
        import cupy as cp

        kernel = cp.RawKernel(
            r"""
            #include <cupy/complex.cuh>
            extern "C" __global__ void primitive(
                const double* u, const complex<double>* c,
                const complex<double>* prefix, complex<double>* out,
                long long n, int panels, int stride) {
                long long i = (long long)blockDim.x*blockIdx.x+threadIdx.x;
                if (i >= n) return;
                double t = (u[i]+.5)*panels;
                int j = min((int)t, panels-1);
                double x = 2*(t-j)-1;
                complex<double> v = c[j*stride+stride-1];
                for (int k=stride-2;k>=0;--k) v=v*x+c[j*stride+k];
                out[i] = prefix[j]+v;
            }""",
            "primitive",
        )
        u = cp.asarray(u, dtype=cp.float64, order="C")
        out = cp.empty(u.shape, dtype=cp.complex128)
        kernel(
            ((u.size + 255) // 256,),
            (256,),
            (
                u,
                cp.asarray(self.coeff, dtype=cp.complex128),
                cp.asarray(self.prefix, dtype=cp.complex128),
                out,
                np.int64(u.size),
                np.int32(self.panels),
                np.int32(self.coeff.shape[1]),
            ),
        )
        return out


class WaveGeometry:
    def __init__(self, q, mode="background", lam=1 / 64, phi=0.37):
        if not (0 <= q < np.pi**2 and 0 < lam < np.inf and np.isfinite(phi)):
            raise ValueError("require 0 <= q < pi^2 and lambda > 0")
        if mode not in {"background", "circular_wave"}:
            raise ValueError("unknown geometry")
        self.q, self.k, self.lam, self.mode = q, np.sqrt(q), lam, mode
        self.factor = (
            2 * lam * np.sqrt(q) / np.sqrt(1 + 4 * q * lam**2) * np.exp(1j * phi)
            if mode == "circular_wave"
            else 0j
        )
        panels = max(256, int(np.ceil(16 / lam)))
        self.primitive = PanelPrimitive(
            lambda u: np.exp(1j * u / lam) / np.cos(self.k * u) ** 2, panels
        )
        self.cdf = PanelPrimitive(self.slice_density, max(2048, panels))
        self.volume = float(np.pi * self.cdf.prefix[-1])

    def m(self, u):
        return np.tan(self.k * np.asarray(u)) / self.k if self.q else np.asarray(u)

    def matrix(self, start, end):
        start, end = np.broadcast_arrays(start, end)
        # The sinc identity avoids cancellation for adjacent events.
        dm = (
            (end - start)
            * np.sinc(self.k * (end - start) / np.pi)
            / (np.cos(self.k * start) * np.cos(self.k * end))
        )
        beta = self.factor * (self.primitive(end) - self.primitive(start))
        return np.stack(
            (dm - beta.real, -beta.imag, -beta.imag, dm + beta.real), axis=-1
        ).reshape(dm.shape + (2, 2))

    def slice_density(self, u):
        left, right = self.matrix(-0.5, u), self.matrix(u, 0.5)
        det = np.linalg.det(left) * np.linalg.det(right) / np.linalg.det(left + right)
        return np.cos(self.k * u) ** 2 * np.sqrt(np.maximum(det, 0))

    def sprinkle(self, n, rng, backend="cpu"):
        """Zero rejection, physical-volume iid sampling; final independent shuffle."""
        if backend == "gpu":
            return self._sprinkle_gpu(n, rng)
        if backend != "cpu":
            raise ValueError("backend must be cpu or gpu")
        target = rng.random(n) * self.cdf.prefix[-1]
        idx = np.searchsorted(self.cdf.prefix, target, side="right") - 1
        lo = idx / self.cdf.panels - 0.5
        hi = lo + self.cdf.width
        fraction = (target - self.cdf.prefix[idx]) / (
            self.cdf.prefix[idx + 1] - self.cdf.prefix[idx]
        )
        u = lo + fraction * self.cdf.width
        for _ in range(6):
            error = self.cdf(u) - target
            hi = np.where(error > 0, u, hi)
            lo = np.where(error < 0, u, lo)
            candidate = u - error / np.maximum(self.slice_density(u), 1e-300)
            u = np.clip(candidate, lo, hi)
        li = np.linalg.inv(self.matrix(-0.5, u))
        ri = np.linalg.inv(self.matrix(u, 0.5))
        # Closed 2x2 SPD square root avoids a large batched eigendecomposition.
        total = li + ri
        d = np.sqrt(np.linalg.det(total))
        root = total + d[:, None, None] * np.eye(2)
        invroot = (
            np.linalg.inv(root)
            * np.sqrt(np.trace(total, axis1=-2, axis2=-1) + 2 * d)[:, None, None]
        )
        s, theta = 1 - np.sqrt(rng.random(n)), rng.uniform(0, 2 * np.pi, n)
        z = np.sqrt(2 * s)[:, None] * np.column_stack((np.cos(theta), np.sin(theta)))
        y = np.einsum("nij,nj->ni", invroot, z)
        bottom = 0.5 * np.einsum("ni,nij,nj->n", y, li, y)
        top = 1 - 0.5 * np.einsum("ni,nij,nj->n", y, ri, y)
        points = np.column_stack((u, rng.uniform(bottom, top), y))
        return points[rng.permutation(n)]

    def _sprinkle_gpu(self, n, rng):
        import cupy as cp

        kernel = cp.RawKernel(
            r"""
        #include <cupy/complex.cuh>
        __device__ complex<double> primitive_at(double u, const complex<double>* c,
            const complex<double>* prefix, int panels, int stride) {
            double t=(u+.5)*panels; int j=min((int)t,panels-1);
            double x=2*(t-j)-1; complex<double> v=c[j*stride+stride-1];
            for(int k=stride-2;k>=0;--k) v=v*x+c[j*stride+k];
            return prefix[j]+v;
        }
        extern "C" __global__ void sprinkle(const double* random, double* out,
            long long n, double q, double factor_re, double factor_im,
            const double* cdf, const double* cc, int nc, int sc,
            const complex<double>* pp, const complex<double>* pc, int np, int sp) {
            long long i=(long long)blockDim.x*blockIdx.x+threadIdx.x;
            if(i>=n) return;
            double target=random[4*i]*cdf[nc]; int lo=0, hi=nc;
            while(hi-lo>1){int mid=(lo+hi)/2;if(cdf[mid]<=target)lo=mid;else hi=mid;}
            int j=lo; double x=2*(target-cdf[j])/(cdf[j+1]-cdf[j])-1;
            double lower=-1, upper=1;
            for(int it=0;it<8;++it){
                double v=cc[j*sc+sc-1], d=0;
                for(int k=sc-2;k>=0;--k){d=d*x+v;v=v*x+cc[j*sc+k];}
                double e=v+cdf[j]-target;
                if(e>0)upper=x; if(e<0)lower=x;
                x=fmin(upper,fmax(lower,x-e/fmax(d,1e-300)));
            }
            double u=-.5+(j+(x+1)/2)/nc, k=sqrt(q);
            double half=q>0?tan(k/2)/k:.5, m=q>0?tan(k*u)/k:u;
            double ml=m+half, mr=half-m;
            complex<double> factor(factor_re,factor_im);
            complex<double> beta=factor*primitive_at(u,pc,pp,np,sp);
            complex<double> end=factor*pp[np], right=end-beta;
            double dl=ml*ml-norm(beta), dr=mr*mr-norm(right);
            double l00=(ml+beta.real())/dl, l01=beta.imag()/dl, l11=(ml-beta.real())/dl;
            double r00=(mr+right.real())/dr, r01=right.imag()/dr;
            double r11=(mr-right.real())/dr;
            double a=l00+r00, b=l01+r01, c=l11+r11, det=sqrt(a*c-b*b);
            double scale=sqrt(a+c+2*det)/((a+det)*(c+det)-b*b);
            double radius=sqrt(2*(1-sqrt(random[4*i+1])));
            double angle=6.2831853071795864769*random[4*i+2];
            double z0=radius*cos(angle), z1=radius*sin(angle);
            double y0=scale*((c+det)*z0-b*z1), y1=scale*((a+det)*z1-b*z0);
            double bottom=.5*(l00*y0*y0+2*l01*y0*y1+l11*y1*y1);
            double top=1-.5*(r00*y0*y0+2*r01*y0*y1+r11*y1*y1);
            out[4*i]=u;out[4*i+1]=bottom+random[4*i+3]*(top-bottom);
            out[4*i+2]=y0;out[4*i+3]=y1;
        }""",
            "sprinkle",
        )
        out = cp.empty((n, 4), dtype=cp.float64)
        kernel(
            ((n + 255) // 256,),
            (256,),
            (
                cp.asarray(rng.random((n, 4))),
                out,
                np.int64(n),
                np.float64(self.q),
                np.float64(self.factor.real),
                np.float64(self.factor.imag),
                cp.asarray(self.cdf.prefix),
                cp.asarray(self.cdf.coeff),
                np.int32(self.cdf.panels),
                np.int32(self.cdf.coeff.shape[1]),
                cp.asarray(self.primitive.prefix),
                cp.asarray(self.primitive.coeff),
                np.int32(self.primitive.panels),
                np.int32(self.primitive.coeff.shape[1]),
            ),
        )
        return cp.asnumpy(out)[rng.permutation(n)]

    def interval_volume(self, first, second, degree=96):
        """Continuum truth oracle; integrates a^2 sqrt(det L det R / det M)."""
        first, second = np.broadcast_arrays(first, second)
        length = second[..., 0] - first[..., 0]
        nodes, weights = legendre.leggauss(degree)
        u = first[..., 0, None] + length[..., None] * (nodes + 1) / 2
        left = self.matrix(first[..., 0, None], u)
        right = self.matrix(u, second[..., 0, None])
        whole = self.matrix(first[..., 0], second[..., 0])
        dy = second[..., 2:] - first[..., 2:]
        margin = (
            second[..., 1]
            - first[..., 1]
            - 0.5 * np.einsum("...i,...ij,...j->...", dy, np.linalg.inv(whole), dy)
        )
        density = np.cos(self.k * u) ** 2 * np.sqrt(
            np.maximum(
                np.linalg.det(left)
                * np.linalg.det(right)
                / np.linalg.det(whole)[..., None],
                0,
            )
        )
        return np.pi * np.maximum(margin, 0) ** 2 * length / 2 * (density @ weights)

    def population(self, points, early, late, degree=96):
        start = np.array([-0.5, 0, 0, 0])
        end = np.array([0.5, 1, 0, 0])
        cols = [
            self.interval_volume(start, points, degree),
            self.interval_volume(points, end, degree),
        ]
        cols.extend(self.interval_volume(a, points, degree) for a in early)
        cols.extend(self.interval_volume(points, b, degree) for b in late)
        return np.stack(cols, axis=-1) / self.volume


class OrderOracle:
    """Only query() and size are handed to the observer; storage stays private.

    Tiles use a fused CUDA kernel (two bool outputs, no float pair temporaries).
    The CPU implementation is an independent vectorized reference.
    """

    def __init__(self, geometry, points, backend="cpu"):
        if backend not in {"cpu", "gpu"}:
            raise ValueError("backend must be cpu or gpu")
        self.size = len(points)
        self.backend = backend
        if backend == "gpu":
            import cupy as cp

            device_points = cp.asarray(points)
            u = device_points[:, 0]
            beta = (
                geometry.factor * geometry.primitive.gpu(u)
                if geometry.factor
                else cp.zeros(len(points), dtype=cp.complex128)
            )
            m = cp.tan(geometry.k * u) / geometry.k if geometry.q else u
            self._data = cp.column_stack((device_points, m, beta.real, beta.imag))
            self._kernel = cp.RawKernel(
                r"""
            extern "C" __global__ void relation(
                const double* d, const long long* a, const long long* b,
                bool* forward, bool* backward, long long na, long long nb) {
                long long i=(long long)blockIdx.x*blockDim.x+threadIdx.x;
                if(i>=na*nb) return;
                const double *x=d+7*a[i/nb], *y=d+7*b[i%nb];
                double du=y[0]-x[0], dv=y[1]-x[1];
                double dx=y[2]-x[2], dy=y[3]-x[3], m=y[4]-x[4];
                double re=y[5]-x[5], im=y[6]-x[6];
                double det=m*m-re*re-im*im;
                double quad=((m+re)*dx*dx+2*im*dx*dy+(m-re)*dy*dy);
                double gap=dv-quad/(2*det);
                forward[i]=(du>0 && gap>=0) ||
                    (du==0 && dx==0 && dy==0 && dv>0);
                backward[i]=(du<0 && gap<=0) ||
                    (du==0 && dx==0 && dy==0 && dv<0);
            }""",
                "relation",
            )
        else:
            beta = geometry.factor * geometry.primitive(points[:, 0])
            self._data = np.column_stack(
                (points, geometry.m(points[:, 0]), beta.real, beta.imag)
            )

    def query(self, rows, columns):
        """Strict relations (rows < columns, columns < rows), shaped (R,X)."""
        rows, columns = (
            np.asarray(rows, dtype=np.int64),
            np.asarray(columns, dtype=np.int64),
        )
        if (
            rows.ndim != 1
            or columns.ndim != 1
            or np.any(rows < 0)
            or np.any(rows >= self.size)
            or np.any(columns < 0)
            or np.any(columns >= self.size)
        ):
            raise ValueError("invalid event indices")
        if rows.size * columns.size > 8192 * 4096:
            raise ValueError("query exceeds tile budget; tile the callback")
        if self.backend == "gpu":
            import cupy as cp

            fw, bw = self._query_gpu(rows, columns)
            return cp.asnumpy(fw), cp.asnumpy(bw)
        x, y = self._data[rows, None], self._data[columns][None]
        du, dv, dx, dy, dm, re, im = np.moveaxis(y - x, -1, 0)
        numerator = (dm + re) * dx**2 + 2 * im * dx * dy + (dm - re) * dy**2
        det = dm**2 - re**2 - im**2
        gap = dv - np.divide(numerator, 2 * det, out=np.zeros_like(dm), where=det != 0)
        equal = (du == 0) & (dx == 0) & (dy == 0)
        return ((du > 0) & (gap >= 0)) | (equal & (dv > 0)), (
            ((du < 0) & (gap <= 0)) | (equal & (dv < 0))
        )

    def _query_gpu(self, rows, columns):
        import cupy as cp

        shape = (len(rows), len(columns))
        fw, bw = cp.empty(shape, dtype=bool), cp.empty(shape, dtype=bool)
        if fw.size:
            self._kernel(
                ((fw.size + 255) // 256,),
                (256,),
                (
                    self._data,
                    cp.asarray(rows, dtype=cp.int64),
                    cp.asarray(columns, dtype=cp.int64),
                    fw,
                    bw,
                    np.int64(len(rows)),
                    np.int64(len(columns)),
                ),
            )
        return fw, bw

    def counts(self, references, probes, early=(), late=()):
        """Optional accelerated callback: order counts ONLY, same split contract."""
        if self.backend != "gpu":
            raise ValueError("accelerated counts require GPU backend")
        import cupy as cp

        groups = [
            np.asarray(g, dtype=np.int64) for g in (references, probes, early, late)
        ]
        if any(g.ndim != 1 or np.any(g < 0) or np.any(g >= self.size) for g in groups):
            raise ValueError("invalid indices")
        references, probes, early, late = groups
        if (
            not len(references)
            or np.intersect1d(references, probes, assume_unique=True).size
        ):
            raise ValueError("references must be nonempty and disjoint from probes")
        result = cp.zeros((len(probes), 2 + len(early) + len(late)), dtype=cp.float64)
        for offset in range(0, len(references), 8192):
            rr = references[offset : offset + 8192]
            # The all-ones channel counts tips in the same exact-integer GEMMs.
            ar = cp.vstack(
                (
                    cp.ones((1, len(rr)), dtype=cp.float32),
                    self._query_gpu(early, rr)[0].astype(cp.float32),
                )
            )
            rb = cp.column_stack(
                (
                    cp.ones(len(rr), dtype=cp.float32),
                    self._query_gpu(rr, late)[0].astype(cp.float32),
                )
            )
            for j in range(0, len(probes), 4096):
                fw, bw = self._query_gpu(rr, probes[j : j + 4096])
                target = result[j : j + fw.shape[1]]
                past = (ar @ fw.astype(cp.float32)).T
                future = bw.astype(cp.float32).T @ rb
                target[:, 0] += past[:, 0]
                target[:, 1] += future[:, 0]
                target[:, 2 : 2 + len(early)] += past[:, 1:]
                target[:, 2 + len(early) :] += future[:, 1:]
        return cp.asnumpy(result) / len(references)
