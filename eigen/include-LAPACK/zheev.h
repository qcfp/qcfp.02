#ifdef __cplusplus
extern "C" void zheev_(
	const char &jobz,		// (input)
	const char &uplo,		// (input)
	const int &n,			// (input)
	complex<double> *a,		// a[n][lda] (input/output)
	const int &lda,			// (input)
	double *w,			// w[n] (output)
	complex<double> *work,		// work[lwork] (workspace/output)
	const int &lwork,		// (input)
	double *rwork,			// rwork[max(1, 3*n-2)] (workspace)
	int &info			// (output)
	);
#else /* ! __cplusplus */
void zheev_(
	const char *jobz,		/* (input) */
	const char *uplo,		/* (input) */
	const int *n,			/* (input) */
	complex_double *a,		/* a[n][lda] (input/output) */
	const int *lda,			/* (input) */
	double *w,			/* w[n] (output) */
	complex_double *work,		/* work[lwork] (workspace/output) */
	const int *lwork,		/* (input) */
	double *rwork,			/* rwork[max(1, 3*n-2)] (workspace) */
	int *info			/* (output) */
	);
#endif /* ! __cplusplus */

