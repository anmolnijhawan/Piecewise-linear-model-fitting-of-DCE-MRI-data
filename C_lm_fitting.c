#include "mex.h"
#include "matrix.h"
#define MA 5
#define NDATA 32
#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}
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


int nl,nh;
float *vector(nl,nh)
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

// Starting the covsrt
void covsrt(covar,ma,lista,mfit)
float **covar;
int ma,lista[],mfit;
{
	int i,j;
	float swap;

	for (j=1;j<ma;j++)
		for (i=j+1;i<=ma;i++) covar[i][j]=0.0;
	for (i=1;i<mfit;i++)
		for (j=i+1;j<=mfit;j++) {
			if (lista[j] > lista[i])
				covar[lista[j]][lista[i]]=covar[i][j];
			else
				covar[lista[i]][lista[j]]=covar[i][j];
		}
	swap=covar[1][1];
	for (j=1;j<=ma;j++) {
		covar[1][j]=covar[j][j];
		covar[j][j]=0.0;
	}
	covar[lista[1]][lista[1]]=swap;
	for (j=2;j<=mfit;j++) covar[lista[j]][lista[j]]=covar[1][j];
	for (j=2;j<=ma;j++)
		for (i=1;i<=j-1;i++) covar[i][j]=covar[j][i];
}



// Starting the gauss function
void gaussj(a,n,b,m)
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
                        nrerror("GAUSSJ: Singular Matrix-1");
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
            nrerror("GAUSSJ: Singular Matrix-2");
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
}


// Starting the Function to be minimized.
void PLfit(float x,float a[],float *y,float dyda[],int ma)
{
    int i;
    *y = 0.0;//printf(" y %f ", *y);float temp;
    for(i=1;i<=5;i++)   //Initializing all the derivatives to zero
        dyda[i] = 0.00;

    if(x<=a[1])
       *y = a[3];
    else if(x>a[1]&&x<=a[2]){
       *y = (a[3]+a[4]*(x-a[1]));
       //printf("%f %f %f %f y =%f",(x-a[1]),(a[4]*(x-a[1])),a[3],(a[3]+a[4]*(x-a[1])),*y);
   //temp= *y;
   }
    else if(x>a[2])
       {*y = a[3]+a[4]*(a[2]-a[1])+a[5]*(x-a[2]);
         //printf("\n \n \n inside a[2] %f\n",*y);
         }
   dyda[3]=1.0;
  // printf("\n out  y= %f",*y);
    if(x>a[1]&& x<=a[2])
    {
//printf("\n in 1 y= %f",*y);

       dyda[1]= -1*a[4];
       dyda[4]=  x-a[1];
  // printf("\n in 2 y= %f",*y);

    }
    else if(x>a[2])
    {//printf("\n in 3 y= %f",*y);

       dyda[1]= -1*a[4];
       dyda[2]= a[4]-a[5];
       dyda[4]= a[2]-a[1];
       dyda[5]= x-a[2];
       //printf("\n in 4  y= %f",*y);

    }

    //printf("inside PLfit temp, x= %f , y=%f \n",x,*y);
    //printf("%f %f %f %f %f \n",dyda[1],dyda[2],dyda[3],dyda[4],dyda[5]);
   // printf("%f %f %f %f %f",a[1],a[2],a[3],a[4],a[5]);

}

//Starting the marcof
void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[],int ma, float **alpha, float beta[], float *chisq,void (*funcs)(float, float [], float *, float [],int ))
{
    int i,j,k,l,m,mfit=0;
    //printf("inside mrqcof \n");
float ymod,wt,sig2i,dy,*dyda;
dyda=vector(1,ma);
for (j=1;j<=ma;j++)
if (ia[j]) mfit++;
for (j=1;j<=mfit;j++) {
        for (k=1;k<=j;k++) alpha[j][k]=0.0;
beta[j]=0.0;
}
//printf("inside mrqcof value of %d %f %f\n",ma,x[1],x[30]);
*chisq=0.0;
for (i=1;i<=NDATA;i++) {
    (*funcs)(x[i],a,&ymod,dyda,ma);
//printf("i= %d,x[i]=%f \n",i,x[i]);
sig2i=1.0/(sig[i]*sig[i]);
dy=y[i]-ymod;

for (j=0,l=1;l<=ma;l++) {
if (ia[l]) {
wt=dyda[l]*sig2i;
for (j++,k=0,m=1;m<=l;m++)
if (ia[m]) alpha[j][++k] += wt*dyda[m];
beta[j] += dy*wt;
}
}

*chisq += dy*dy*sig2i;   //find x squared
//
//printf(" inside marcof chisq=%f",*chisq);
}
for (j=2;j<=mfit;j++)     //Fill in the symmetric side.
for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
  free_vector(dyda,1,ma);
}




// Starting the mrqmin function
void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],int ma, float **covar, float **alpha, float *chisq,void (*funcs)(float, float [], float *, float [],int ), float *alamda,int counti)
{

float tempa;
int i,j,k,l;
static int mfit;
static float ochisq,*atry,*beta,*da,**oneda;
if (*alamda < 0.0) {
      atry=vector(1,ma);
beta=vector(1,ma);
da=vector(1,ma);
for (mfit=0,j=1;j<=ma;j++)
if (ia[j]) mfit++;

oneda=matrix(1,mfit,1,1);
*alamda=0.001;
//printf("going to mrqcof\n");
mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
//printf("Coming back from first mrqcof call\n");
//printf("a[1]=%f, a[2]=%f \n",a[1],a[2]);
//printf("chisq %f,ochisq %f\n",*chisq,ochisq);
//printf("Printing alpha \n");
ochisq=(*chisq);
//printf("printing the parameters a[] \n");
//for(i=1;i<=5;i++)
  //  printf("a[%d]=%f\n",i,a[i]);
for (j=1;j<=ma;j++) atry[j]=a[j];
}
for (j=1;j<=mfit;j++) { //Alter linearized fitting matrix, by augmenting difor
for(k=1;k<=mfit;k++) covar[j][k]=alpha[j][k]; //agonal elements.
covar[j][j]=alpha[j][j]*(1.0+(*alamda));
oneda[j][1]=beta[j];
}
//printf("chisq =%f\n",*chisq);

gaussj(covar,mfit,oneda,1);   //Matrix solution.
for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
//printf("printing da[]\n");
//for(i=1;i<=5;i++)
   // printf("da[%d]=%f\n",i,da[i]);

if (*alamda == 0.0) { //Once converged, evaluate covariance matrix.
covsrt(covar,ma,ia,mfit);
covsrt(alpha,ma,ia,mfit); //Spread out alpha to its full size too.
free_matrix(oneda,1,mfit,1,1);
free_vector(da,1,ma);
free_vector(beta,1,ma);
free_vector(atry,1,ma);
return;
}
for (j=0,l=1;l<=ma;l++) //Did the trial succeed?
if (ia[l]) atry[l]=a[l]+da[++j];
//mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
//printf("Checking the Success\n");
//printf("atry[1]=%f, atry[2]=%f",atry[1],atry[2]);
//printf("Chisq =%f Ochisq=%f \n",*chisq,ochisq);
//printf("atry[1]=%f, atry[2]=%f\n",atry[1],atry[2]);


//&&(atry[1]>=5.0)&&(atry[1]<=9.0)&&(atry[2]>=7)&&(atry[2]<12.0)&&(atry[1]<atry[2])

ochisq = (*chisq);
if((atry[1]<5.0)||(atry[1]>9.0)||(atry[2]<7.0)||(atry[2]>12.0)||(atry[1]>atry[2]))
    *chisq = ochisq;
else
{
    mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
}
if ((*chisq < ochisq)) { //Success, accept the new solution.
//printf("\n Success \n");
*alamda *= 0.1;

//ochisq=(*chisq);
for (j=1;j<=mfit;j++) {
for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
beta[j]=da[j];
}
for (l=1;l<=ma;l++) a[l]=atry[l];

if((a[1]<5.0)||(a[1]>9.0)||(a[1]>a[2]))
    a[1] = 5.5;
if((a[2]<7.0)||(a[2]>12.0)||(a[1]>a[2]))
    a[2] = 8.5;
if((a[2]-a[1])<=1.0)
     a[1] = a[1]-1.0;
if(a[3]<0.0||a[3]>=0.05)
        a[3]= 0.0005;
if(a[4]<=0.0||a[4]>1.0)
            a[4]=0.45;
if(a[5]<=-1.0||a[5]>=1.0)
      a[5]=-0.005;

}
else {    //Failure, increase alamda and return.
if(counti>=4)
{
 //  printf("\n");
   //printf("\n");
    *alamda *= 15.0;
    if((a[1]<5.0)||(a[1]>9.0)||(a[1]>a[2]))
    a[1] = 5.5;
if((a[2]<7.0)||(a[2]>12.0)||(a[1]>a[2]))
    a[2] = 8.5;

if((a[2]-a[1])<=1.0)
     a[1] = a[1]-1.0;
   if(a[3]<0.0||a[3]>=0.05)
     a[3]=0.0005;
   if(a[4]<=0.0||a[4]>1.0)
            a[4]=0.45;
if(a[5]<=-1||a[5]>=1)
      a[5]=-0.005;

}
 else{
*alamda *= 10.0;
if((a[1]<5.0)||(a[1]>9.0)||(a[1]>a[2]))
    a[1] = 5.5;
if((a[2]<7.0)||(a[2]>12.0)||(a[1]>a[2]))
    a[2] = 8.5;

if((a[2]-a[1])<1.0)
     a[1] = a[1]-1.0;

   if(a[3]<0.0||a[3]>+0.05)
     a[3]=0.0005;
   if(a[4]<=0.0||a[4]>1.0)
            a[4]=0.45;
if(a[5]<=-1.0||a[5]>=1.0)
      a[5]=-0.005;

 }
//printf("printing the value of alamda=%f\n",*alamda);
//printf("printing the values of chisq=%f ,ochisq=%f \n",*chisq,ochisq);
//*chisq=ochisq;   // checkthis line..

//printf("New values\n");
//printf("a[1]=%f a[2]\n%f",a[1],a[2]);
}
//printf("Printing the alpha matrix");

//for(i=1;i<=MA;i++)
   // {
       // for(j=1;j<=MA;j++)
        //printf("%f ",alpha[i][j]);
        //printf("\n");
    //}
//printf("%p \n",funcs);
//tempa= *chisq;
}



// Start of program.
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
   int M;double *v;double *ss;int zz;
    int count=2;float new_a[6];

  //printf("Hello \n");
  //printf("Hello!!!\n");
    float *x,*y,*sig,**covar,**alpha,*a,chisq,ochisq,alambda;
    int *ia,ma;
    char name[80],initname[80],outname[80];
    //FILE *in,*init,*out;
    int i,k,ndata,j,jj;
    x= vector(1,NDATA);

    y= vector(1,NDATA);

    sig= vector(1,NDATA);
    //if(sig==NULL)
      //  printf("Null \n");
    a = vector(1,MA);
    covar= matrix(1,MA,1,MA);
    alpha=matrix(1,MA,1,MA);
    ia=ivector(1,MA);
    ma=MA;

    M = mxGetM(prhs[0]);
    v= mxGetPr(prhs[0]);
    for(i=1;i<=NDATA;i++)
 {
    y[i]=v[i-1];
 }

    for(i=1;i<=NDATA;i++)
    {
        x[i]= 0.065*i;
       sig[i]=1.00;
    }


a[1]=(4.0*0.065);  //Initializing the parameters
  a[2]=(7.0*0.065);   // a[1]=alpha, a[2]=beta, a[3]=c(constant)
  a[3]=0.0;  // a[4]=b1(slope 1)  a[5] = b2(slope 2)
  a[4]=0.5;
  a[5]=0.00;   //All five parameters initialized

    //for(i=1;i<=MA;i++){
      //  for(j=1;j<=MA;j++){
        //   covar[i][j]=i+j;
         //  alpha[i][j]=i+j;}}
 for(i=1;i<=MA;i++)
    {
        ia[i]=1;
    }
 // printf("printing the alambda \n");
  alambda = -1.11;
   chisq = 0.0;
   //printf("calling mrqmin \n");
     mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,&chisq,PLfit,&alambda,1);
     printf("chisq =%f\n",chisq);
     for(jj=1;jj<=MA;jj++)    // Initializing the jacobian matrix with zeroes.
{
    for(k=1;k<=MA;k++)
        printf("%f ",covar[jj][k]);
        printf("\n");
}
  printf("\n");

    //  printf("Printing the values after mrqmin \n");
      //printf(" PLfit %p\n",PLfit);
  //printf("alambda %f \n",alambda);
 //printf("chisq %f \n", chisq);
 //ochisq = chisq + 2.0;
 //while(fabs((ochisq - chisq)/ochisq) > 0.000001 ) {
//printf("1 a[1]=%f a[2]=%f a[3]=%f a[4]=%f a[5]=%f chisq=%f\n",(a[1]*0.065),(a[2]*0.065),a[3],a[4],a[5],chisq);
   //printf("chisq %f \n", chisq);
/*mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,&chisq,PLfit,&alambda,2);
//printf("2 a[1]=%f a[2]=%f a[3]=%f a[4]=%f a[5]=%f  chisq=%f\n",(a[1]*0.065),(a[2]*0.065),a[3],a[4],a[5],chisq);
 // printf("chisq %f \n", chisq);
     ochisq=chisq;
     new_a[1]=a[1];
         new_a[2]=a[2];
         new_a[3]=a[3];
         new_a[4]=a[4];
         new_a[5]=a[5];

 while(count<=30)
 {
   mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,&chisq,PLfit,&alambda,(count-1));
//printf("3 a[1]=%f a[2]=%f a[3]=%f a[4]=%f a[5]=%f chisq=%f\n",(a[1]*0.065),(a[2]*0.065),a[3],a[4],a[5],chisq);

     if(ochisq > (chisq))
     {

         new_a[1]=a[1];
         new_a[2]=a[2];
         new_a[3]=a[3];
         new_a[4]=a[4];
         new_a[5]=a[5];
         ochisq = chisq ;
     }
     count++;
 }*/

//printf("4 a[1]=%f a[2]=%f a[3]=%f a[4]=%f a[5]=%f chisq=%f\n",(new_a[1]*0.065),(new_a[2]*0.065),new_a[3],new_a[4],new_a[5],ochisq);
//printf("\n%f %f\n",new_a[1],new_a[2]);
//FILE *f3;
/*f3 = fopen("best_fit.txt","w");
   for(i=1;i<=32;i++)
			{
				if(i<=new_a[1])
				{
				  fprintf(f3,"%f\n",new_a[3]);
				}
				else if(i>new_a[1]&&i<=new_a[2])
				{
					fprintf(f3,"%f\n",(new_a[3]+new_a[4]*(i-new_a[1])));
				}
				else if(i>new_a[2])
				{
					fprintf(f3,"%f\n",(new_a[3]+(new_a[4]*(new_a[2]-new_a[1]))+(new_a[5]*(i-new_a[2]))));

				}
			}
			fclose(f3);*/

  /* for(i=1;i<=MA;i++)
    {
        //for(j=1;j<=MA;j++)
        printf("%d ",ia[i]);
        printf("%f  ",a[i]);
        printf("\n");
    }
// printf("%f ",sig[2]);
    for(i=1;i<=NDATA;i++)
    {
        printf("%f ",x[i]);
        printf("%f ",sig[i]);
        printf("%f ",y[i]);
        printf("\n");
    }
    for(i=1;i<=MA;i++)
    {
        for(j=1;j<=MA;j++)
        printf("%f ",covar[i][j]);
        printf("\n");
    }
    for(i=1;i<=MA;i++)
    {
        for(j=1;j<=MA;j++)
        printf("%f ",alpha[i][j]);
        printf("\n");
    }*/
   plhs[0]= mxCreateDoubleMatrix(6, 1, mxREAL);
     ss= mxGetPr(plhs[0]);
     ss[0]=new_a[1];
	 ss[1]=new_a[2];
	 ss[3]=new_a[4];
	 ss[4]=new_a[5];




    free_vector(x,1,NDATA);
    free_vector(y,1,NDATA);
    free_vector(a,1,MA);
    free_vector(sig,1,NDATA);
    free_matrix(covar,1,MA,1,MA);
    free_matrix(alpha,1,MA,1,MA);
    free_vector(ia,1,MA);


}
