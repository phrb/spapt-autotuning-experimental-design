/*@ begin PerfTuning (
  def build
  {
    arg build_command = 'timeout --kill-after=30s --signal=9 20m gcc -O3 -fopenmp -DDYNAMIC';
    arg libs = '-lm -lrt';
  }

  def performance_counter
  {
    arg repetitions = 10;
  }

  def performance_params
  {
    param T_I[] = [1,2,4,8,16,32,64,128,256,512,1024,2048];
    param T_J[] = [1,2,4,8,16,32,64,128,256,512,1024,2048];

    param RT_I[] = [1,2,4,8,16,32];
    param RT_J[] = [1,2,4,8,16,32];

    param U_I[] = range(1,31);
    param U_J[] = range(1,31);

    param SCREP[] = [False,True];
    param VEC[] = [False,True];
    # param OMP[] = [False,True];

    constraint reg_capacity = (4*U_I*U_J + 3*U_I + 3*U_J <= 128);
    constraint unroll_limit = ((U_I == 1) or (U_J == 1));
  }

  def search
  {
    arg algorithm = 'GPR';
    arg total_runs = 1;
    arg gpr_starting_sample = 13;
    arg gpr_steps = 31;
    arg gpr_extra_experiments = 12;
    arg gpr_testing_set_size = 300000;
    arg gpr_failure_multiplier = 80;
  }

  def input_params
  {
    let RANGE = 8000;
    let BSIZE = 512*32;
    param SIZE = RANGE;
    param N = RANGE;
  }

  def input_vars
  {
    decl static double X0[N][N] = random;
    decl static double X1[N][N] = random;
    decl static double X2[N][N] = random;
    decl static double Y[N][N] = 0;
    decl static double u0[N] = random;
    decl static double u1[N] = random;
    decl static double u2[N] = random;
    decl double a0 = 32.12;
    decl double a1 = 3322.12;
    decl double a2 = 1.123;
    decl double b00 = 1321.9;
    decl double b01 = 21.55;
    decl double b02 = 10.3;
    decl double b11 = 1210.313;
    decl double b12 = 9.373;
    decl double b22 = 1992.31221;
  }

  def validation {
    arg validation_file = 'validation_3x.c';
  }
) @*/

int i,j,ii,jj,iii,jjj,it,jt;

#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))


/*@ begin Loop(
  transform Composite(
    tile = [('i',T_I,'ii'),('j',T_J,'jj')],
    unrolljam = (['i','j'],[U_I,U_J]),
    scalarreplace = (SCREP, 'double', 'scv_'),
    regtile = (['i','j'],[RT_I,RT_J]),
    vector = (VEC, ['ivdep','vector always']),
    openmp = (False, 'omp parallel for private(i,j,ii,jj,iii,jjj,Y,X0,X1,X2)')
  )

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    {
      Y[i][j]=a0*X0[i][j] + a1*X1[i][j] + a2*X2[i][j]
    + 2.0*b00*u0[i]*u0[j]
    + 2.0*b11*u1[i]*u1[j]
    + 2.0*b22*u2[i]*u2[j]
    + b01*u0[i]*u1[j] + b01*u1[i]*u0[j]
    + b02*u0[i]*u2[j] + b02*u2[i]*u0[j]
    + b12*u1[i]*u2[j] + b12*u2[i]*u1[j];
    }

) @*/

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    {
      Y[i][j]=a0*X0[i][j] + a1*X1[i][j] + a2*X2[i][j]
    + 2.0*b00*u0[i]*u0[j]
    + 2.0*b11*u1[i]*u1[j]
    + 2.0*b22*u2[i]*u2[j]
    + b01*(u0[i]*u1[j] + u1[i]*u0[j])
    + b02*(u0[i]*u2[j] + u2[i]*u0[j])
    + b12*(u1[i]*u2[j] + u2[i]*u1[j]);
    }

/*@ end @*/
/*@ end @*/
