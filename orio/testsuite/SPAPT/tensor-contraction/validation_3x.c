
int isValid() {

  double actual = 878906250000.000000;

  double s_sum = 0.0;
  double q_sum = 0.0;
  double rand1=0.1, rand2=0.9;
  double expected=0.0;
  int i,j,k,l,m;
  double diff=0.0;


  for(i=0; i<=V-1; i++)
    for(j=0; j<=V-1; j++)
      for(k=0; k<=O-1; k++)
        for(l=0; l<=O-1; l++)
      for(m=0; m<=O-1; m++)
        s_sum+=R[i][j][k][l];


  expected = s_sum;

  diff=abs(expected-actual);

  fprintf(stderr, "expected=%f\n",expected);
  fprintf(stderr, "actual=%f\n",actual);
  fprintf(stderr, "diff=%f\n",diff);
  fprintf(stderr, "diff=%d\n",(diff < 0.00000001));

  if (diff < 0.00000001)
    return 1;
  else
    return 0;
}




