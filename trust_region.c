#include "mex.h"
#include "matrix.h"
#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}

//--------------------------------------------------------------------
//Functions for memory allocation and deallocation
//--------------------------------------------------------------------
void nrerror(error_text)
char error_text[];
{

	printf("%s\n",error_text);
}
void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
	free((char*) (v+nl));
}

float *vector(int nl,int nh)
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}
float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	float **m;

	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}
void free_ivector(v,nl,nh)
int *v,nl,nh;
{
	free((char*) (v+nl));
}

int *ivector(nl,nh)
int nl,nh;
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}
//------------------------------------------------------------------------

float signum(float a)
{
    if(a<0.0)
        return -1.0;
    else if(a>0.0)
        return 1.0;
    else if(a==0.0)
        return 0.0;
}

float maximum(float a, float b)
{
    if(a>b)
        return a;
    if(b>a)
        return b;
    if(a==b)
        return a;
}

float minimum(float a, float b)
{
    if(a<b)
        return a;
    if(b<a)
        return b;
    if(b==a)
        return a;
}

//-----------------------------------------------------------------------------
//Gauss Jordan elimination
//-----------------------------------------------------------------------------
int gaussj(a,n,b,m)
float **a,**b;
int n,m;
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll,*ivector();
	float big,dum,pivinv;
	void nrerror(),free_ivector();

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++)
        ipiv[j]=0;
	for (i=1;i<=n;i++)
    {
		big=0.0;
		for (j=1;j<=n;j++)
        {
			if (ipiv[j] != 1)
            {
				for (k=1;k<=n;k++)
                {
					if (ipiv[k] == 0)
                    {
						if (fabs(a[j][k]) >= big)
                        {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
                    else if (ipiv[k] > 1)
                    {
                      nrerror("GAUSSJ: Singular Matrix-1");
                      return 1;
                    }
                        				}
            }
        }
		++(ipiv[icol]);
		if (irow != icol)
        {
			for (l=1;l<=n;l++)
                SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++)
                SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0)
        {
            nrerror("GAUSSJ: Singular Matrix-2");
            return 1;
        }
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++)
            a[icol][l] *= pivinv;
		for (l=1;l<=m;l++)
            b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
        {
			if (ll != icol)
            {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
        }
	}
	for (l=n;l>=1;l--)
    {
		if (indxr[l] != indxc[l])
        {
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
	return 0;
}


//-------------------------------------------------------------------------------------------------
//Function for calculating the values at all time points according to Piece-wise Linear Function
//-------------------------------------------------------------------------------------------------
void getF(float x[],float a[],int NDATA,float F[])
{
   int i;
   for(i=1;i<=NDATA;i++)
   {
    if(x[i]<=a[1])
      F[i]=a[3];
    else if((x[i]>a[1])&&(x[i]<=a[2]))
      F[i]=(a[3]+(a[4]*(x[i]-a[1])));
    else if(x[i]>a[2])
      F[i]=(a[3]+(a[4]*(a[2]-a[1]))+(a[5]*(x[i]-a[2])));
   }
}

void chisq(float x[],float y[],float a[],int NDATA,float *chisq1)
{
    int i;float *F;
    F=vector(1,NDATA);
    getF(x,a,NDATA,F);
    *chisq1=0.0;
    for(i=1;i<=NDATA;i++)
    {
       *chisq1+=(F[i]-y[i])*(F[i]-y[i]);
    }
    free_vector(F,1,NDATA);
}

//-----------------------------------------------------------------------------
// Get the vector v assigned according to gradients.
//-----------------------------------------------------------------------------
void getV(float dF[],float a[],float v[],float u[],float l[],int MA)
{
   int i;
   for(i=1;i<=MA;i++)
   {
       if(dF[i]<0.0)
       {
           v[i]=a[i]-u[i];
       }
       else if(dF[i]>=0.0)
       {
           v[i]=a[i]-l[i];
       }
   }
}

//--------------------------------------------------------------------------------------
//Function for calculating the gradient of the minimizing function and Hessian matrix.
//--------------------------------------------------------------------------------------
void gradient(float x[],float y[],float a[],int NDATA,int MA,float dF[],float **hessian)
{
    int i,j,k;float *f;float **deltaF;
    deltaF=matrix(1,MA,1,NDATA);
    f=vector(1,NDATA);
    getF(x,a,NDATA,f);
    for(i=1;i<=MA;i++)
    {
        for(j=1;j<=NDATA;j++)
        {
            deltaF[i][j]=0.0;
        }
        dF[i]=0.0;
    }
    for(i=1;i<=MA;i++)
    {
        for(j=1;j<=MA;j++)
        {
            hessian[i][j]=0.0;
        }
    }

    for(i=1;i<=NDATA;i++)
    {
        if(x[i]<=a[1])
        {
            deltaF[3][i]=1.0;
        }
        else if((x[i]>a[1])&&(x[i]<=a[2]))
        {
            deltaF[3][i]=1.0;
            deltaF[1][i]=(-1*a[4]);
            deltaF[4][i]=(x[i]-a[1]);
        }
        else if(x[i]>a[2])
        {
            deltaF[3][i]=1.0;
            deltaF[1][i]=(-1*a[4]);
            deltaF[2][i]=(a[4]-a[5]);
            deltaF[4][i]=(a[2]-a[1]);
            deltaF[5][i]=(x[i]-a[2]);
        }
    }
    for(j=1;j<=MA;j++)
    {
      for(i=1;i<=NDATA;i++)
      {
        dF[j]+=((f[i]-y[i])*deltaF[j][i]);
      }
    }
    for(i=1;i<=MA;i++)
    {
        for(j=1;j<=MA;j++)
        {
            for(k=1;k<=NDATA;k++)
            {
                hessian[i][j]+=deltaF[i][k]*deltaF[j][k];
            }
        }
    }
    free_matrix(deltaF,1,MA,1,NDATA);
    free_vector(f,1,NDATA);
}


void getINVERSEDk(float v[],float d[],int MA)
{
    int i;
    for(i=1;i<=MA;i++)
    {
        d[i]=sqrt(fabs(v[i]));
    }
}

//----------------------------------------------------------------------------
//Getting value of determinant to see if matrix is positive semi-definite
//-----------------------------------------------------------------------------
void getCofactor(float **mat,float **temp,int p,int q,int n)
{
    int i,j,row,col;
     i=1; j=1;

    // Looping for each element of the matrix
    for (row=1;row<=n;row++)
    {
        for (col=1;col<=n;col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = mat[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == n )
                {
                    j = 1;
                    i++;
                }
            }
        }
    }
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of mat[][]. */
float determinantOfMatrix(float **mat, int n)
{
    float D,sign;float **temp;int f;
    temp=matrix(1,5,1,5);
     D = 0.0; // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return mat[1][1];

    //int temp[N][N]; // To store cofactors
    //printf("printing mat[][]\n");
    /*for(i=1;i<=n;i++)
    {
        for(j=1;j<=n;j++)
        {
            printf("%f ",mat[i][j]);
        }printf("\n");
    }*/
    sign = 1.0;  // To store sign multiplier

     // Iterate for each element of first row
    for (f = 1; f<=n; f++)
    {
        // Getting Cofactor of mat[0][f]
        getCofactor(mat, temp, 1, f, n);
        D += sign * mat[1][f] * determinantOfMatrix(temp, n-1);

        // terms are to be added with alternate sign
        sign = -sign;
    }
    free_matrix(temp,1,5,1,5);
    return D;

}
//----------------------------------------------------------------------------------

float check_positive_definite(int MA,float **M)
{
    int i;float temp,flag;
    flag=0.0;
    for(i=1;i<=MA;i++)
    {
          temp=determinantOfMatrix(M, i);
          //printf("%d determinant= %f\n",i,temp);
          if(temp<=0.0)
          {
              flag=1.0;  //if determinant is negative then flag set to 1
              break;
          }
    }
    return flag;
}
//----------------------------------------------------------------------------------

//---------------------------------------------------------------------------
//method for finding norm
//------------------------------------------------------------------
float get_norm(float n[],int M)
{
    int i;float temp;
    temp=0.0;
    for(i=1;i<=M;i++)
    {
        temp+=n[i]*n[i];
    }
     temp=sqrt(temp);
     return temp;
}

//-------------------------------------------------------------------------------------
//Cholesky decomposition for reducing the matrix into lower and upper triangular matix
//-------------------------------------------------------------------------------------

int choldc(float **a, int n, float p[],float **r)
{
      void nrerror(char error_text[]);
      int i,j,k;
      float sum;
      for (i=1;i<=n;i++) {
         for (j=i;j<=n;j++) {
           for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
           if (i == j) {
              if (sum < 0.0){   //a, with rounding errors, is not positive definite.
                  nrerror("choldc failed");
                  return 1;
                  }
              p[i]=sqrt(sum);
              r[i][i]=sqrt(sum);
           } else
           {
             a[j][i]=sum/p[i];
             r[j][i]=sum/p[i];
           }
         }
       }
    return 0;
}
//--------------------------------------------------------------------------------------
//Function for solving Ax=b based on cholesky decomposition.
//------------------------------------------------------------------------------------------
void cholsl(float **a, int n, float p[], float b[], float x[])
{
  int i,k;
  float sum;
  for (i=1;i<=n;i++)
  {                  //Solve L.y = b, storing y in x.
     for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
     x[i]=sum/p[i];
  }
  for (i=n;i>=1;i--)
  {                     //Solve LT.x = y.
     for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
     x[i]=sum/p[i];
  }
}


int solving_step_2(float **B,float g[],int MA,float delta,float HAT_s[])
{
    int count;float check,lambda,new_lambda,norm_p_L,norm_q_L,**tempB,**R_T,**Rt,*p_L,*q_L,*p,*x,**p_L1,norm,norm1;int i,j,l,k;float flag;
    tempB=matrix(1,MA,1,MA);
    p=vector(1,MA);
    x=vector(1,MA);
    R_T=matrix(1,MA,1,MA);
    Rt=matrix(1,MA,1,MA);
    p_L=vector(1,MA);
    q_L=vector(1,MA);
    p_L1=matrix(1,MA,1,1);
    lambda=0.0;
    flag=0.0;

    for(l=1;l<=MA;l++)
    {
       for(k=1;k<=MA;k++)
       {
          tempB[l][k]=B[l][k];
          R_T[l][k]=0.0;
       }
    }
    flag=check_positive_definite(MA,tempB);
    //printf("flag= %f\n",flag);
    if(flag==0.0)
    {
        //printf("First matrix is semi-positive\n");
        choldc(tempB,MA,p,R_T);
        cholsl(tempB,MA,p,g,x);
        norm=get_norm(x,MA);

     //printf("lambda=%f and norm=%f and delta=%f \n",lambda,norm,delta);
      if((lambda==0.0)&&(norm<=delta))
      {
         //printf("Entering if\n");
        for(l=1;l<=MA;l++)
            HAT_s[l]=(x[l]*(-1));
      } //printf("HAT_s[]\n");  for(l=1;l<=MA;l++) printf("%f ",HAT_s[l]);
      else if(norm>delta)
      {
        //printf("norm > delta \n");

         count=1;check=1;
         while((check>=0.02))
         {
             count++;
             if(count>12)
                return 2;
             /*for(l=1;l<=MA;l++)
             {
               for(k=1;k<=MA;k++)
               {
                 tempB[l][k]=B[l][k];  //Re-initializing the matrix because cholesky decomposition destroyed the matrix.
               }
             }
             for(l=1;l<=MA;l++)
             {
              tempB[l][l]=(tempB[l][l]*(1+lambda));
             }*/
              i=choldc(tempB,MA,p,R_T);  //put some condition to see if it can be decomposed into lower triangular matrix or not
              if(i==1)
                return 1;
              cholsl(tempB,MA,p,g,x);
              norm1=get_norm(x,MA);
           //printf("printing lower triangular matrix\n");
           /*for(l=1;l<=MA;l++)
           {
             for(k=1;k<=MA;k++)
             {
                 printf("%f ",R_T[l][k]);
             }printf("\n");
           }*/
           //printf("norm1=%f delta=%f \n",norm1,delta);
           for(i=1;i<=MA;i++)
           {
               p_L[i]=(x[i]*(-1));
               p_L1[i][1]=(x[i]*(-1));
           }
           norm_p_L=get_norm(p_L,MA);  //Get norm of p_L
           i=gaussj(R_T,MA,p_L1,1);  //for getting q_L
           if(i==1)
                return 1;
              //printf("Completing gauss\n");
           for(i=1;i<=MA;i++)
             q_L[i]=p_L1[i][1];
          norm_q_L=get_norm(q_L,MA);
          //printf("norm_p=%f norm_q=%f delta=%f \n",norm_p_L,norm_q_L,delta);
          new_lambda=lambda+(((norm_p_L/norm_q_L)*(norm_p_L/norm_q_L))*((norm_p_L-delta)/delta));
          //printf("new_lambda=%f, lambda=%f delta=%f\n",new_lambda,lambda,delta);
          lambda=new_lambda;

          for(l=1;l<=MA;l++)
             {
               for(k=1;k<=MA;k++)
               {
                 tempB[l][k]=B[l][k];  //Re-initializing the matrix because cholesky decomposition destroyed the matrix.
               }
             }
             for(l=1;l<=MA;l++)
             {
              tempB[l][l]=(tempB[l][l]*(1+lambda));
             }

             flag=check_positive_definite(MA,tempB);
             //printf("flag=%f\n",flag);
             if(flag==0.0)
             {
                 for(l=1;l<=MA;l++)
                   HAT_s[l]=(x[l]*(-1));
                 //printf("Count =%d\n",count);
             }
             //printf("lambda=%f delta=%f norm_p_L=%f\n",lambda,delta,norm_p_L);
             check=(norm_p_L-delta)/delta;
         }

      }
    }
    //printf("getting out\n");
    free_matrix(tempB,1,MA,1,MA);
    free_vector(p,1,MA);
    free_vector(x,1,MA);
    free_matrix(R_T,1,MA,1,MA);
    free_matrix(Rt,1,MA,1,MA);
    free_vector(p_L,1,MA);
    free_vector(q_L,1,MA);
    free_matrix(p_L1,1,MA,1,1);
    return 0;
}
//updating_trust_region(rho,mu,delta_k,&delta_k1);
void updating_trust_region(float rho,float mu,float delta_k,float *delta_k1)
{
    float eta,delta_l,gamma1,gamma2,gamma_O,lambda_l,lambda_u;
    eta=0.75;  gamma1=0.5; gamma2=2.0; gamma_O=0.0625;
    lambda_u=sqrt(56);
    lambda_u=maximum(lambda_u,1.0);
    //printf("lambda_u=%f\n",lambda_u);
    if(rho<=0.0)
    {
        *delta_k1=gamma_O*delta_k;
    }
    if((rho>0.0)&&(rho<=mu))
    {
        *delta_k1=(gamma1*delta_k)/2.0;
    }
    if((rho>mu)&&(rho<eta))
    {
       *delta_k1=(delta_k-((delta_k-(gamma1*delta_k))/2.0));
    }
    if(rho>=eta)
    {
        if(delta_k>delta_l)
        {
            *delta_k1=((delta_k+(gamma2*delta_k))/2.0);
        }
        else
            *delta_k1=((delta_k+(gamma2*delta_k))/2.0);
    }


}

void updating_trust_region1(float rho,float mu,float delta_k,float *delta_k1)
{
    float eta,delta_l,gamma1,gamma2,gamma_O,lambda_l,lambda_u;
    eta=0.75;  gamma1=0.5; gamma2=2.0; gamma_O=0.0625; delta_l=1.0;
    lambda_u=sqrt(56);
    lambda_u=maximum(lambda_u,1.0);
    //printf("lambda_u=%f\n",lambda_u);
    if(rho<=0.0)
    {
        *delta_k1=gamma_O*delta_k;
    }
    if((rho>0.0)&&(rho<=mu))
    {
        *delta_k1=(gamma1*delta_k);
    }
    if((rho>mu)&&(rho<eta))
    {
       *delta_k1=(delta_k);
    }
    if(rho>=eta)
    {
        if(delta_k>delta_l)
        {
            *delta_k1=(gamma2*delta_k);

        }
        else
            *delta_k1=minimum((gamma2*delta_k),lambda_u);
            //printf("(gamma2*delta_k)= %f  lambda_u=%f \n",(gamma2*delta_k),lambda_u);
            //printf("delta k1=%f\n",(*delta_k1));
    }


}


void main_iteration(float x[],float y[],float a[],int NDATA,int MA,float u[],float l[])
{
    float mu,psi,*a_n,*new_x;int i,j,flag1;float *a_new,delta,delta_k,delta_k1,**hessian,*dF,*Jv,*v,*d,error,chisq1,chisq2,*s,temp,temp1,*HAT_s,rho,**HAT_M,*HAT_g;
    int count;
    a_new=vector(1,MA);
    a_n=vector(1,MA);
    s=vector(1,MA);
    HAT_s=vector(1,MA);
    HAT_M=matrix(1,MA,1,MA);
    hessian=matrix(1,MA,1,MA);    // hessian matrix
    dF=vector(1,MA);   // gradient vector
    Jv=vector(1,MA);
    HAT_g=vector(1,MA);
    v=vector(1,MA);          // vector v for building matrix D
    d=vector(1,MA);          // Diagonal matrix D
    count=1;
    gradient(x,y,a,NDATA,MA,dF,hessian);
    //delta=0.3;
    delta_k=minimum((0.1*get_norm(dF,MA)),7.4833);
    //printf("initial delta=%f\n",delta);

    while(count<=20)
   {
       count++;
   //step 1
    chisq(x,y,a,NDATA,&error); //printf("error from initial guess %f\n",*error);  //computing f_k
    gradient(x,y,a,NDATA,MA,dF,hessian);  // computing gradient and hessian
    /*printf("printing hessian matrix[][]\n");
    for(i=1;i<=MA;i++)
    {
       for(j=1;j<=MA;j++)
       {
         printf("%f ",hessian[i][j]);
       } printf("\n");
    }*/
    /*printf("printing dF\n");
    for(i=1;i<=MA;i++)
      printf("%f ",dF[i]);
    printf("\n");
    printf("norm of g=%f\n",get_norm(dF,MA));*/

    getV(dF,a,v,u,l,MA);  // Computing vector v depending on gradient
    /*printf("printing vector v\n");
    for(i=1;i<=MA;i++)
     printf("%f ",v[i]);
    printf("\n");*/

    getINVERSEDk(v,d,MA);  // Computing inverse of matrix D
    /*printf("printing vector d\n");
    for(i=1;i<=MA;i++)
     printf("%f ",d[i]);
    printf("\n");*/

    //step 2 with consideration that (Xk + Sk) belong to int(F)
     //transforming the problem to HAT parameters.Original minimization problem transformed into new minimization problem.
     for(i=1;i<=MA;i++)
    {
       Jv[i]=signum(dF[i]);   //Getting the matrix Jv[][].
    }
    /*printf("printing the Jv matrix\n");
    for(i=1;i<=MA;i++)
     printf("%f ",Jv[i]);
    printf("\n");*/

    for(i=1;i<=MA;i++)
    {
         HAT_g[i]=d[i]*dF[i];  //filling the vector HAT_g[].
    }
    /*printf("printing HAT_g\n");
    for(i=1;i<=MA;i++)
     printf("%f ",HAT_g[i]);
    printf("\n");*/

    for(i=1;i<=MA;i++)
    {
        for(j=1;j<=MA;j++)    //building of matrix HAT_M[][].
        {
            if(i==j)
            HAT_M[i][j]=(d[i]*hessian[i][j]*d[j])+(dF[i]*Jv[i]);
            else
            HAT_M[i][j]=(d[i]*hessian[i][j]*d[j]);
        }
    }
    /*printf("printing matrix HAT_M[][] \n");
    for(i=1;i<=MA;i++)
    {
        for(j=1;j<=MA;j++)    //building of matrix HAT_M[][].
        {
           printf("%f ",HAT_M[i][j]);
        } printf("\n");
    }*/

     //step 2
     //printf("delta=%f\n",delta_k);
     //printf("getting in of solving_step_2\n");
     i=solving_step_2(HAT_M,HAT_g,MA,delta_k,HAT_s);
     //printf("getting out of solving_step_2\n");
     //printf("i=%d",i);
     if((i==0)||(i==2))
     {
        for(i=1;i<=MA;i++)
     {
        s[i]=d[i]*HAT_s[i];   // Conversion from HAT_s to s by multiplying HAT_s with matrix D inverse.
     }
     for(i=1;i<=MA;i++)
     {
       //printf("i=%d %f \n",i,(a[i]+s[i]));
       a_n[i]=(a[i]+s[i]);   //new parameters are assigned accordingly.
     }
     flag1=0;
     //printf("Checking the boundaries\n");
     for(i=1;i<=MA;i++)
     {
        if((a_n[i]>=l[i])&&(a_n[i]<=u[i]))
        {
            flag1++;
           // printf("%f ",a_n[i]);
        }
        //printf("\n");
     }
     //printf("flag=%d\n",flag1);

    if(flag1==5)
    {
      //step 3
    for(i=1;i<=MA;i++)
    {
        s[i]=d[i]*HAT_s[i];   // Conversion from HAT_s to s by multiplying HAT_s with matrix D inverse.
    }
    //printf("\n in main s[1]=%f s[2]=%f s[3]=%f s[4]=%f s[5]=%f \n",s[1],s[2],s[3],s[4],s[5]);
    for(i=1;i<=MA;i++)
    {
       //printf("i=%d %f \n",i,(a[i]+s[i]));
       a_n[i]=(a[i]+s[i]);   //new parameters are assigned accordingly.
    }
    chisq(x,y,a_n,NDATA,&error);
    //printf("new_error=%f\n",error);
    /*printf("new a[]\n");
    for(i=1;i<=MA;i++)
    {
      printf("%f ",a_n[i]);
    }printf("\n");*/

    //building of psi_HAT
    temp=0.0;
    temp1=0.0;
    for(i=1;i<=MA;i++)
    {
        temp+=(HAT_g[i]*HAT_s[i]);
        for(j=1;j<=MA;j++)
        {
            temp1+=HAT_s[i]*HAT_M[i][j]*HAT_s[j];
        }
    }
     psi=temp+(temp1/2.0);    //psi is the value of quadratic model.
     chisq1=0.0;
    chisq(x,y,a_n,NDATA,&chisq1);
    //printf("chisq1=%f \n",chisq1);
    chisq(x,y,a,NDATA,&chisq2);
    temp=0.0;
    for(i=1;i<=MA;i++)
    {
      temp+=HAT_s[i]*dF[i]*Jv[i]*HAT_s[i];
    }
    temp=temp/2.0;
    rho=(chisq1-chisq2+temp)/(psi);  // value of quadratic expression is still missing
    //printf("rho= %f chisq1= %f chisq2= %f psi= %f\n",rho,chisq1,chisq2,psi);

    //step 4
    mu=0.25;
    if(rho>mu)
    {
        //printf("inside if step 4\n");
       for(i=1;i<=MA;i++)
       a[i]=a[i]+s[i];  // not very sure about s[i] check once.
    }
    else
    {
        //printf("Inside else step 4\n");
        for(i=1;i<=MA;i++)
        a[i]=a[i];
    }
    /*printf("checking updated parametrs\n");
    for(i=1;i<=MA;i++)
        printf("%f ",a[i]);
    printf("\n");*/

    //step 5
    //printf("delta_k=%f rho=%f\n",delta_k,rho);
    updating_trust_region(rho,mu,delta_k,&delta_k1);
     //printf("delta_k1=%.7f \n",delta_k1);
     delta_k=delta_k1;
     temp=(a[2]-a[1])/0.065;
     if(temp>=3.7)
     {
        a[1]=a[1]+0.065;
        a[2]=a[2]-0.065;
     }

    }
   }
   if(i==1)
   {
     a[1]=7.2;
     a[2]=9.8;
     delta_k=1.10*delta_k;
   }
   // go for singular matrix ones.
 }

 // printf("flag=%d\n",flag1);
    free_vector(a_new,1,MA);
    free_vector(a_n,1,MA);
    free_vector(dF,1,MA);
    free_matrix(hessian,1,MA,1,MA);
    free_matrix(HAT_M,1,MA,1,MA);
    free_vector(s,1,MA);
    free_vector(v,1,MA);
    free_vector(d,1,MA);
    free_vector(HAT_s,1,MA);
    free_vector(Jv,1,MA);
    free_vector(HAT_g,1,MA);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 int M;double *v;double *ss;int NDATA,i,MA;float chisq11,new_chisq,onset,*y,*a,*x,*u,*l;
 M = mxGetM(prhs[0]);
 v= mxGetPr(prhs[0]);
 NDATA = (M);   // NDATA represents the no. of concentration values given as input to mex function
 MA=5;          // MA represents the number of parameters to optimized.
 y= vector(1,NDATA);  // Allocating space for storing Concentration values.
 a = vector(1,MA);    // Allocating space for storing parameter values.
 x = vector(1,NDATA);  // Allocating space for storing time point values.
 u = vector(1,MA);    // for storing the upper bound for parameters to be optimized
 l = vector(1,MA);    // for storing the lower bound for parameters to be optimized
 for(i=1;i<=NDATA;i++)   // For extracting the data and putting it in y[] array.
  {
     y[i]=v[i-1];              // For storing concentration values
     x[i]=(i-1)*0.065;        // For storing the time points.
      //printf("%f ",x[i]);
      //printf("%f \n",y[i]);
  } //printf("\n");
 onset =0.0;            // Initializing the constant parameter
 for(i=1;i<=6;i++)
   onset += y[i];
 onset = onset/6.0;    // Initializing the constant parameter with average of first six concentration values.
 a[1]=(6.5*0.065);     // Array a[] represents the five parameters to be fit, they are BAT, Beta, constant, slope1 and slope2
 a[2]=(9.5*0.065);
 a[3]=onset;              // Initializing the parameters
 //printf("a[3]=%f\n",a[3]);
 a[4]=0.1;             // These values represent the initial guesses.
 a[5]=0.00;
 u[1]=(9.0*0.065); u[2]=(12.0*0.065); u[3]=(onset+(0.5*onset)); u[4]=1.0; u[5]=1.0;  //Assigning the upper-bound for parameter values.
 l[1]=(4.5*0.065); l[2]=(7.0*0.065); l[3]=(onset-(0.5*onset)); l[4]=0.0; l[5]=-1.0;  //Assigning the lower-bound for parameter values.
 chisq11 = 0.00;     //
 new_chisq = 0.0;
 chisq(x,y,a,NDATA,&chisq11);
 //printf("chisq=%f \n",chisq11);
 main_iteration(x,y,a,NDATA,MA,u,l);

 plhs[0]= mxCreateDoubleMatrix(5, 1, mxREAL);  //Assigning space for output parameters from mex-file.
 ss= mxGetPr(plhs[0]);  //getting the pointer of output parameter from mex-file
 ss[0]=a[1];      // Optimized parameter values are returned to matlab by mex files
 ss[1]=a[2];
 ss[2]=a[3];
 ss[3]=a[4];
 ss[4]=a[5];


 free_vector(y,1,NDATA);
 free_vector(x,1,NDATA);
 free_vector(a,1,MA);
 free_vector(u,1,MA);
 free_vector(l,1,MA);

}

