/*
 * hmcox.c -- pathwise coordinate descent for the hierarchical penalized
 * stratified Cox model used in IPD meta-analysis.
 *
 *   minimise  Q(alpha, eps~) = L(eta)
 *                              + sum_j P(|alpha_j|; lambda_a)
 *                              + sum_{k,j} P(|eps~_kj|; lambda_e)
 *   subject to  alpha_j = 0  =>  eps~_kj = 0 for all k   (hierarchy)
 *
 *   L(eta) = -(1/N) * stratified log partial likelihood (Breslow ties),
 *   eta_i  = sum_j x_ij * (alpha_j + eps~_{k(i) j} / omega_k),
 *   omega_k = sqrt(n_k / N)  (standardisation of the study-specific columns).
 *
 * Rows of X must be sorted by study and, within study, by follow-up time.
 * Each outer iteration builds the diagonal quadratic surrogate of L at the
 * current eta (h_i = haz_i * A_i / N, Simon et al. 2011; Breheny & Huang 2011),
 * minimises the penalised surrogate by exact block coordinate descent that
 * never leaves the hierarchy set, and accepts the step only if the true
 * penalised objective does not increase (step-halving otherwise).
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* penalties                                                          */
/* ------------------------------------------------------------------ */
static double pen_val(double b, double lam, double gam, int type) {
  double a = fabs(b);
  if (lam <= 0 || a == 0) return 0.0;
  if (type == 0) return lam * a;                                 /* lasso */
  if (type == 1)                                                 /* MCP   */
    return (a <= gam * lam) ? lam * a - a * a / (2.0 * gam) : 0.5 * gam * lam * lam;
  /* SCAD */
  if (a <= lam) return lam * a;
  if (a <= gam * lam) return (2.0 * gam * lam * a - a * a - lam * lam) / (2.0 * (gam - 1.0));
  return lam * lam * (gam + 1.0) / 2.0;
}

static double sgn(double x) { return (x > 0) - (x < 0); }

/* exact global minimiser of  f(b) = v/2 (b - z)^2 + P(|b|; lam, gam) */
static double pen_min(double z, double v, double lam, double gam, int type) {
  if (lam <= 0) return z;
  double cand[8];
  int nc = 0;
  double s = sgn(z), vz = v * fabs(z);
  cand[nc++] = 0.0;
  cand[nc++] = z;
  if (type == 0) {
    cand[nc++] = (vz > lam) ? s * (vz - lam) / v : 0.0;
  } else if (type == 1) {
    if (v > 1.0 / gam) {
      double b = (vz > lam) ? s * (vz - lam) / (v - 1.0 / gam) : 0.0;
      if (fabs(b) <= gam * lam) cand[nc++] = b;
    }
    cand[nc++] = s * gam * lam;
  } else {
    double b1 = (vz > lam) ? s * (vz - lam) / v : 0.0;
    if (fabs(b1) <= lam) cand[nc++] = b1;
    if (v > 1.0 / (gam - 1.0)) {
      double t = gam * lam / (gam - 1.0);
      double b2 = (vz > t) ? s * (vz - t) / (v - 1.0 / (gam - 1.0)) : 0.0;
      if (fabs(b2) > lam && fabs(b2) <= gam * lam) cand[nc++] = b2;
    }
    cand[nc++] = s * lam;
    cand[nc++] = s * gam * lam;
  }
  double best = 0.0, fbest = R_PosInf;
  for (int c = 0; c < nc; c++) {
    double b = cand[c];
    double f = 0.5 * v * (b - z) * (b - z) + pen_val(b, lam, gam, type);
    if (f < fbest - 1e-15) { fbest = f; best = b; }
  }
  return best;
}

/* ------------------------------------------------------------------ */
/* stratified partial likelihood pieces                               */
/* ------------------------------------------------------------------ */
typedef struct {
  int N, p, K;
  const double *X, *d, *omega;
  const int *sstart, *tfirst, *tlast, *sid;
  double *haz, *R, *C, *rsk, *A, *buf;
} coxdat;

/* negative log partial likelihood / N ; optionally fills h and r */
static double cox_irls(coxdat *D, const double *eta, double *h, double *r) {
  double loss = 0.0;
  int N = D->N;
  for (int k = 0; k < D->K; k++) {
    int a = D->sstart[k], b = D->sstart[k + 1];
    if (b <= a) continue;
    double m = eta[a];
    for (int i = a; i < b; i++) if (eta[i] > m) m = eta[i];
    for (int i = a; i < b; i++) D->haz[i] = exp(eta[i] - m);
    D->R[b - 1] = D->haz[b - 1];
    for (int i = b - 2; i >= a; i--) D->R[i] = D->R[i + 1] + D->haz[i];
    for (int i = a; i < b; i++) D->rsk[i] = D->R[D->tfirst[i]];
    double c = 0.0;
    for (int i = a; i < b; i++) {
      if (D->d[i] > 0) c += D->d[i] / D->rsk[i];
      D->C[i] = c;
    }
    for (int i = a; i < b; i++) {
      D->A[i] = D->C[D->tlast[i]];
      if (D->d[i] > 0) loss -= D->d[i] * (eta[i] - m - log(D->rsk[i]));
    }
    if (h != NULL) {
      for (int i = a; i < b; i++) {
        double ha = D->haz[i] * D->A[i];
        h[i] = ha / N;
        r[i] = (ha > 0) ? (D->d[i] - ha) / ha : 0.0;
      }
    }
  }
  return loss / N;
}

static void compute_eta(coxdat *D, const double *a, const double *e, double *eta) {
  int N = D->N, K = D->K;
  for (int i = 0; i < N; i++) eta[i] = 0.0;
  for (int j = 0; j < D->p; j++) {
    const double *x = D->X + (size_t)j * N;
    if (a[j] != 0) for (int i = 0; i < N; i++) eta[i] += x[i] * a[j];
    for (int k = 0; k < K; k++) {
      double ek = e[k + (size_t)K * j];
      if (ek == 0) continue;
      double c = ek / D->omega[k];
      for (int i = D->sstart[k]; i < D->sstart[k + 1]; i++) eta[i] += x[i] * c;
    }
  }
}

static double pen_total(int p, int K, const double *a, const double *e,
                        double lamA, double lamE, double gam, int type) {
  double s = 0.0;
  for (int j = 0; j < p; j++) {
    s += pen_val(a[j], lamA, gam, type);
    for (int k = 0; k < K; k++) s += pen_val(e[k + (size_t)K * j], lamE, gam, type);
  }
  return s;
}

/* one hierarchy-feasible block update for covariate j; returns max change */
static double block_update(coxdat *D, int j, double *a, double *e, double *h,
                           double *r, double *eta, double lamA, double lamE,
                           double gam, int type, int het) {
  int N = D->N, K = D->K;
  const double *x = D->X + (size_t)j * N;
  double maxch = 0.0;
  double xwr = 0.0, xwx = 0.0;
  for (int i = 0; i < N; i++) { double hx = h[i] * x[i]; xwr += hx * r[i]; xwx += hx * x[i]; }
  if (xwx <= 1e-14) return 0.0;
  double v = xwx;
  double anew = pen_min(xwr / v + a[j], v, lamA, gam, type);

  int epsnz = 0;
  if (het) for (int k = 0; k < K; k++) if (e[k + (size_t)K * j] != 0) { epsnz = 1; break; }

  if (anew == 0.0 && epsnz) {
    /* infeasible to keep deviations with a zero global effect:
       compare the current block with the all-zero block */
    double dq = 0.0, dp = -pen_val(a[j], lamA, gam, type);
    for (int k = 0; k < K; k++) {
      double ek = e[k + (size_t)K * j];
      dp -= pen_val(ek, lamE, gam, type);
      double c = a[j] + ek / D->omega[k];
      for (int i = D->sstart[k]; i < D->sstart[k + 1]; i++) {
        double de = -x[i] * c;
        dq += h[i] * (-r[i] * de + 0.5 * de * de);
      }
    }
    if (dq + dp < 0) {
      for (int k = 0; k < K; k++) {
        double c = a[j] + e[k + (size_t)K * j] / D->omega[k];
        for (int i = D->sstart[k]; i < D->sstart[k + 1]; i++) {
          double de = -x[i] * c;
          r[i] -= de; eta[i] += de;
        }
        double ch = fabs(e[k + (size_t)K * j]);
        if (ch > maxch) maxch = ch;
        e[k + (size_t)K * j] = 0.0;
      }
      if (fabs(a[j]) * sqrt(v) > maxch) maxch = fabs(a[j]) * sqrt(v);
      a[j] = 0.0;
    }
    return maxch;
  }

  double shift = anew - a[j];
  if (shift != 0) {
    for (int i = 0; i < N; i++) { double si = shift * x[i]; r[i] -= si; eta[i] += si; }
    a[j] = anew;
    maxch = fabs(shift) * sqrt(v);
  }
  if (het && a[j] != 0) {
    for (int k = 0; k < K; k++) {
      /* sharing rule as a hard constraint: the last study that shares
         alpha_j keeps a zero deviation */
      if (e[k + (size_t)K * j] == 0) {
        int nzero = 0;
        for (int kk = 0; kk < K; kk++) if (e[kk + (size_t)K * j] == 0) nzero++;
        if (nzero <= 1) continue;
      }
      int s0 = D->sstart[k], s1 = D->sstart[k + 1];
      double om = D->omega[k];
      double wr = 0.0, wx = 0.0;
      for (int i = s0; i < s1; i++) { double hx = h[i] * x[i]; wr += hx * r[i]; wx += hx * x[i]; }
      wr /= om; wx /= (om * om);
      if (wx <= 1e-14) continue;
      size_t idx = k + (size_t)K * j;
      double enew = pen_min(wr / wx + e[idx], wx, lamE, gam, type);
      double sh = enew - e[idx];
      if (sh != 0) {
        double c = sh / om;
        for (int i = s0; i < s1; i++) { double si = c * x[i]; r[i] -= si; eta[i] += si; }
        e[idx] = enew;
        if (fabs(sh) * sqrt(wx) > maxch) maxch = fabs(sh) * sqrt(wx);
      }
    }
    /* Re-centering move: (alpha_j + c, eps_kj - c) leaves every theta_kj and
       hence the likelihood unchanged; choose c to minimise the penalty. The
       penalty is concave (MCP/SCAD) or linear (lasso) in c between the
       breakpoints c = eps_kj, so the minimum is attained at a breakpoint.
       c = -alpha_j is excluded because it would violate the hierarchy.
       If all K deviations are nonzero, the best breakpoint is taken even if it
       does not lower the penalty: the identification rule requires at least
       one study to share alpha_j (otherwise alpha_j and eps_.j are not
       separately identified and the refit design is singular). */
    int nnz = 0;
    for (int k = 0; k < K; k++) if (e[k + (size_t)K * j] != 0) nnz++;
    if (nnz > 0) {
      double best_c = 0.0, best_p = pen_val(a[j], lamA, gam, type);
      double *orig = D->buf;
      int KK = K;
      for (int k = 0; k < KK; k++) {
        orig[k] = e[k + (size_t)K * j] / D->omega[k];      /* raw scale */
        best_p += pen_val(e[k + (size_t)K * j], lamE, gam, type);
      }
      for (int c0 = 0; c0 < KK; c0++) {
        double c = orig[c0];
        if (c == 0 || a[j] + c == 0) continue;
        double pv = pen_val(a[j] + c, lamA, gam, type);
        for (int k = 0; k < KK; k++) pv += pen_val(D->omega[k] * (orig[k] - c), lamE, gam, type);
        if (pv < best_p - 1e-14 || (nnz == K && best_c == 0)) { best_p = pv; best_c = c; }  /* forced only if infeasible */
      }
      if (best_c != 0) {
        a[j] += best_c;
        for (int k = 0; k < K; k++) {
          e[k + (size_t)K * j] = D->omega[k] * (orig[k] - best_c);
        }
        for (int k = 0; k < K; k++) if (orig[k] == best_c) e[k + (size_t)K * j] = 0.0;
        if (fabs(best_c) * sqrt(v) > maxch) maxch = fabs(best_c) * sqrt(v);
      }
    }
  }
  return maxch;
}

static void project_hier(int p, int K, const double *a, double *e) {
  for (int j = 0; j < p; j++)
    if (a[j] == 0) for (int k = 0; k < K; k++) e[k + (size_t)K * j] = 0.0;
}

/* projection onto the sharing rule: for any covariate whose K deviations are
   all nonzero, move to the likelihood-preserving representation with the
   smallest penalty among those with one zero deviation */
static void project_share(coxdat *D, double *a, double *e, double lamA, double lamE,
                          double gam, int type) {
  int K = D->K;
  for (int j = 0; j < D->p; j++) {
    if (a[j] == 0) continue;
    int nnz = 0;
    for (int k = 0; k < K; k++) if (e[k + (size_t)K * j] != 0) nnz++;
    if (nnz < K) continue;
    double best_p = R_PosInf, best_c = 0.0;
    for (int c0 = 0; c0 < K; c0++) {
      double c = e[c0 + (size_t)K * j] / D->omega[c0];
      if (a[j] + c == 0) continue;
      double pv = pen_val(a[j] + c, lamA, gam, type);
      for (int k = 0; k < K; k++)
        pv += pen_val(D->omega[k] * (e[k + (size_t)K * j] / D->omega[k] - c), lamE, gam, type);
      if (pv < best_p) { best_p = pv; best_c = c; }
    }
    if (best_p < R_PosInf) {
      double *orig = D->buf;
      for (int k = 0; k < K; k++) orig[k] = e[k + (size_t)K * j] / D->omega[k];
      a[j] += best_c;
      for (int k = 0; k < K; k++)
        e[k + (size_t)K * j] = (orig[k] == best_c) ? 0.0 : D->omega[k] * (orig[k] - best_c);
    }
  }
}

/* ------------------------------------------------------------------ */
/* main path routine                                                  */
/* ------------------------------------------------------------------ */
SEXP hmc_path(SEXP X_, SEXP d_, SEXP sstart_, SEXP tfirst_, SEXP tlast_,
              SEXP omega_, SEXP type_, SEXP gamma_, SEXP lam_, SEXP ratio_,
              SEXP het_, SEXP tol_, SEXP maxouter_, SEXP maxinner_,
              SEXP dfmax_, SEXP a0_, SEXP e0_) {
  int N = length(d_);
  int p = length(X_) / N;
  int K = length(omega_);
  int L = length(lam_);
  int type = INTEGER(type_)[0], het = INTEGER(het_)[0];
  double gam = REAL(gamma_)[0], ratio = REAL(ratio_)[0], tol = REAL(tol_)[0];
  int max_outer = INTEGER(maxouter_)[0], max_inner = INTEGER(maxinner_)[0];
  int dfmax = INTEGER(dfmax_)[0];
  double *lam = REAL(lam_);

  coxdat D;
  D.N = N; D.p = p; D.K = K;
  D.X = REAL(X_); D.d = REAL(d_); D.omega = REAL(omega_);
  D.sstart = INTEGER(sstart_); D.tfirst = INTEGER(tfirst_); D.tlast = INTEGER(tlast_);
  D.haz = (double *)R_alloc(N, sizeof(double));
  D.R = (double *)R_alloc(N, sizeof(double));
  D.C = (double *)R_alloc(N, sizeof(double));
  D.rsk = (double *)R_alloc(N, sizeof(double));
  D.A = (double *)R_alloc(N, sizeof(double));
  D.buf = (double *)R_alloc(K, sizeof(double));

  double *a = (double *)R_alloc(p, sizeof(double));
  double *e = (double *)R_alloc((size_t)K * p, sizeof(double));
  double *a_old = (double *)R_alloc(p, sizeof(double));
  double *e_old = (double *)R_alloc((size_t)K * p, sizeof(double));
  double *a_new = (double *)R_alloc(p, sizeof(double));
  double *e_new = (double *)R_alloc((size_t)K * p, sizeof(double));
  double *eta = (double *)R_alloc(N, sizeof(double));
  double *h = (double *)R_alloc(N, sizeof(double));
  double *r = (double *)R_alloc(N, sizeof(double));
  int *act = (int *)R_alloc(p, sizeof(int));
  for (int j = 0; j < p; j++) a[j] = REAL(a0_)[j];
  for (size_t q = 0; q < (size_t)K * p; q++) e[q] = het ? REAL(e0_)[q] : 0.0;
  project_hier(p, K, a, e);
  compute_eta(&D, a, e, eta);

  SEXP A_out = PROTECT(allocMatrix(REALSXP, p, L));
  SEXP E_out = PROTECT(allocVector(REALSXP, (size_t)K * p * L));
  SEXP loss_out = PROTECT(allocVector(REALSXP, L));
  SEXP obj_out = PROTECT(allocVector(REALSXP, L));
  SEXP iter_out = PROTECT(allocVector(INTSXP, L));
  SEXP conv_out = PROTECT(allocVector(INTSXP, L));
  SEXP half_out = PROTECT(allocVector(INTSXP, L));
  double *Ao = REAL(A_out), *Eo = REAL(E_out);
  for (size_t q = 0; q < (size_t)p * L; q++) Ao[q] = NA_REAL;
  for (size_t q = 0; q < (size_t)K * p * L; q++) Eo[q] = NA_REAL;
  for (int l = 0; l < L; l++) {
    REAL(loss_out)[l] = NA_REAL; REAL(obj_out)[l] = NA_REAL;
    INTEGER(iter_out)[l] = NA_INTEGER; INTEGER(conv_out)[l] = 0; INTEGER(half_out)[l] = 0;
  }

  for (int l = 0; l < L; l++) {
    R_CheckUserInterrupt();
    double lamA = lam[l], lamE = ratio * lam[l];
    int conv = 0, nhalf = 0, outer;
    double loss = 0.0;
    for (outer = 1; outer <= max_outer; outer++) {
      memcpy(a_old, a, p * sizeof(double));
      memcpy(e_old, e, (size_t)K * p * sizeof(double));
      double loss_old = cox_irls(&D, eta, h, r);
      double obj_old = loss_old + pen_total(p, K, a, e, lamA, lamE, gam, type);

      /* inner: exact block coordinate descent on the surrogate with active sets */
      int it = 0;
      while (it < max_inner) {
        double mc = 0.0;
        for (int j = 0; j < p; j++) {
          double c = block_update(&D, j, a, e, h, r, eta, lamA, lamE, gam, type, het);
          if (c > mc) mc = c;
        }
        it++;
        for (int j = 0; j < p; j++) act[j] = (a[j] != 0);
        if (mc < tol) break;
        while (it < max_inner) {
          double mc2 = 0.0;
          for (int j = 0; j < p; j++) {
            if (!act[j]) continue;
            double c = block_update(&D, j, a, e, h, r, eta, lamA, lamE, gam, type, het);
            if (c > mc2) mc2 = c;
          }
          it++;
          if (mc2 < tol) break;
        }
      }

      /* true objective at the candidate; step-halving if it went up */
      loss = cox_irls(&D, eta, NULL, NULL);
      double obj_new = loss + pen_total(p, K, a, e, lamA, lamE, gam, type);
      if (obj_new > obj_old + 1e-10 * (1.0 + fabs(obj_old))) {
        memcpy(a_new, a, p * sizeof(double));
        memcpy(e_new, e, (size_t)K * p * sizeof(double));
        double t = 0.5;
        int accepted = 0;
        for (int hh = 0; hh < 30; hh++) {
          for (int j = 0; j < p; j++) a[j] = a_old[j] + t * (a_new[j] - a_old[j]);
          for (size_t q = 0; q < (size_t)K * p; q++) e[q] = e_old[q] + t * (e_new[q] - e_old[q]);
          project_hier(p, K, a, e);
          if (het) project_share(&D, a, e, lamA, lamE, gam, type);
          compute_eta(&D, a, e, eta);
          loss = cox_irls(&D, eta, NULL, NULL);
          double obj_t = loss + pen_total(p, K, a, e, lamA, lamE, gam, type);
          nhalf++;
          if (obj_t <= obj_old + 1e-10 * (1.0 + fabs(obj_old))) { accepted = 1; break; }
          t *= 0.5;
        }
        if (!accepted) {                 /* no ascent direction found: stop */
          memcpy(a, a_old, p * sizeof(double));
          memcpy(e, e_old, (size_t)K * p * sizeof(double));
          compute_eta(&D, a, e, eta);
          loss = cox_irls(&D, eta, NULL, NULL);
          conv = 1;
          break;
        }
      }
      double md = 0.0;
      for (int j = 0; j < p; j++) if (fabs(a[j] - a_old[j]) > md) md = fabs(a[j] - a_old[j]);
      for (size_t q = 0; q < (size_t)K * p; q++) if (fabs(e[q] - e_old[q]) > md) md = fabs(e[q] - e_old[q]);
      if (md < tol) { conv = 1; break; }
    }
    /* store (deviations returned on the original scale) */
    for (int j = 0; j < p; j++) Ao[(size_t)l * p + j] = a[j];
    for (int j = 0; j < p; j++)
      for (int k = 0; k < K; k++)
        Eo[(size_t)l * K * p + k + (size_t)K * j] = e[k + (size_t)K * j] / D.omega[k];
    REAL(loss_out)[l] = loss;
    REAL(obj_out)[l] = loss + pen_total(p, K, a, e, lamA, lamE, gam, type);
    INTEGER(iter_out)[l] = (outer > max_outer) ? max_outer : outer;
    INTEGER(conv_out)[l] = conv;
    INTEGER(half_out)[l] = nhalf;
    int nact = 0;
    for (int j = 0; j < p; j++) if (a[j] != 0) nact++;
    if (nact > dfmax) break;
  }

  SEXP res = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(res, 0, A_out);
  SET_VECTOR_ELT(res, 1, E_out);
  SET_VECTOR_ELT(res, 2, loss_out);
  SET_VECTOR_ELT(res, 3, obj_out);
  SET_VECTOR_ELT(res, 4, iter_out);
  SET_VECTOR_ELT(res, 5, conv_out);
  SET_VECTOR_ELT(res, 6, half_out);
  UNPROTECT(8);
  return res;
}

/* stratified negative log partial likelihood / N and score X'(d - haz*A)/N */
SEXP hmc_loss_score(SEXP X_, SEXP d_, SEXP sstart_, SEXP tfirst_, SEXP tlast_,
                    SEXP eta_, SEXP want_score_) {
  int N = length(d_);
  int p = (length(X_) > 0) ? length(X_) / N : 0;
  int K = length(sstart_) - 1;
  coxdat D;
  D.N = N; D.p = p; D.K = K;
  D.X = (p > 0) ? REAL(X_) : NULL; D.d = REAL(d_);
  D.sstart = INTEGER(sstart_); D.tfirst = INTEGER(tfirst_); D.tlast = INTEGER(tlast_);
  D.haz = (double *)R_alloc(N, sizeof(double));
  D.R = (double *)R_alloc(N, sizeof(double));
  D.C = (double *)R_alloc(N, sizeof(double));
  D.rsk = (double *)R_alloc(N, sizeof(double));
  D.A = (double *)R_alloc(N, sizeof(double));
  double loss = cox_irls(&D, REAL(eta_), NULL, NULL);
  int ws = INTEGER(want_score_)[0];
  SEXP res = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, ScalarReal(loss));
  if (ws && p > 0) {
    SEXP sc = PROTECT(allocVector(REALSXP, p));
    double *s = REAL(sc);
    for (int j = 0; j < p; j++) {
      const double *x = D.X + (size_t)j * N;
      double acc = 0.0;
      for (int i = 0; i < N; i++) acc += x[i] * (D.d[i] - D.haz[i] * D.A[i]);
      s[j] = acc / N;
    }
    SET_VECTOR_ELT(res, 1, sc);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return res;
}

static const R_CallMethodDef CallEntries[] = {
  {"hmc_path", (DL_FUNC)&hmc_path, 17},
  {"hmc_loss_score", (DL_FUNC)&hmc_loss_score, 7},
  {NULL, NULL, 0}};

void R_init_hiermetacox(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
