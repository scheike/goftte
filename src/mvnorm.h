template <matrix_order O1, matrix_style S1,
  matrix_order O2, matrix_style S2,
  matrix_order O3, matrix_style S3>
  double dmvnorm (const Matrix<double, O1, S1>& x,
		  const Matrix<double, O2, S2>& mu,
		  const Matrix<double, O3, S3>& Sigma)
{
  return ( std::exp(lndmvn(x,mu,Sigma)) );
}


