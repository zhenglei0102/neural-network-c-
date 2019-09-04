//Matrix process

#ifndef _MATRIX_PROCESS_
#define _MATRIX_PROCESS_

#include "afxwin.h"
#include "strstrea.h"
#include "fstream.h"
#include "iomanip.h"
#include "stdlib.h"
#include "malloc.h"
#include "math.h"

template<class Type>
class CMatrix :public CObject
{
	public:		
			Type    **f;
			int     row;
			int     col;
	public:
		int  GetRow(){return row;}	//Row numbers
		int  GetCol(){return col;}	//colum numbers
		void SetRow(int r){row=r;}
		void SetCol(int l){col=l;}
		Type** GetF(){return f;}	
		CMatrix();
		CMatrix(int,int);
		CMatrix(CMatrix &copy);
		CMatrix(Type *input, int r, int c);
		CMatrix(Type **input, int r, int c);
		~CMatrix();

		BOOL Ones(int);
		BOOL Ones(int,int);
		BOOL Zeros(int);
		BOOL Zeros(int,int);

		CMatrix Turn();	
		CMatrix Inverse();

		CMatrix Construct(CMatrix &mat,CString style);
		CMatrix GetChild(int,int,int,int);
		CMatrix Outer(CMatrix x,CMatrix y);
		BOOL Diag(Type *,int);
		BOOL Malloc(int,int);
		BOOL IsInverse();
    
		CMatrix   operator*(CMatrix  &m);
		friend CMatrix  operator*(double,CMatrix);
		friend CMatrix  operator*(CMatrix,double);
		CMatrix   operator/(Type x);
		CMatrix   operator+(CMatrix  &m);
		CMatrix   operator+(Type x);
		CMatrix   operator-(CMatrix  &m);
		CMatrix   operator-(Type x);
		CMatrix   operator-();
		CMatrix   operator~(){return Turn();}
		CMatrix   operator!(){return Inverse();}
		Type* &operator[](int);
		friend istream  &operator>>(istream  &,CMatrix  &);
		friend ostream  &operator<<(ostream  &,CMatrix  &);
		void operator=(const CMatrix  &m);
		void operator=(Type);
	
		Type ** MatrixAlloc(int,int);	//old
		void    MatrixFree(Type **);    //old

		void fMatrixAlloc(int r,int c); //new
		void fMatrixFree();				//new

		void    Hessenberg();
		BOOL    Cholesky(); 
		int rank();	
		CMatrix QREigen(double,int);	
		CMatrix Jacobi(double eps,int jt);
		//int Svd(Type *u, Type *v, double eps);
};

template<class Type>
Type* &CMatrix<Type>::operator[](int i)
{
	Type *tp;
	tp=(Type *)f[i];
	return tp;
}
template<class Type>
CMatrix<Type>::CMatrix()
{
	f=NULL;
	row=0;
	col=0;
}
template<class Type>
CMatrix<Type>::CMatrix(int m,int n)
{	
	fMatrixAlloc(m,n);
	row=m;
	col=n;  
}
template<class Type>
CMatrix<Type>::CMatrix(Type **input, int r, int c)
{
	fMatrixAlloc(r,c);
	row=r;
	col=c;
	for(int i=0;i<r;i++)
	for(int j=0;j<c;j++)
	{
		f[i][j]=input[i][j];
	}
}

template<class Type>
CMatrix<Type>::CMatrix(Type *input, int r, int c)
{
	fMatrixAlloc(r,c);
	row=r;
	col=c;
	for(int i=0;i<r;i++)
	for(int j=0;j<c;j++)
	{
		f[i][j]=input[i*c+j];
	}
}
template<class Type>
CMatrix<Type>::~CMatrix()
{
	if(f!=NULL) fMatrixFree();
//	if(f!=NULL) MatrixFree(f);
}

template<class Type>
BOOL CMatrix<Type>::Ones(int m)
{ 	
	fMatrixAlloc(m,m);
	for(int i=0;i<m;i++)
	for(int j=0;j<m;j++)
	{
		if(i!=j)f[i][j]=0;
		else f[i][j]=1;
	}
	row=m;
	col=m;
	return TRUE;
}
template<class Type>
BOOL CMatrix<Type>::Ones(int m,int n)
{ 	
	fMatrixAlloc(m,n);
	for(int i=0;i<m;i++)
	for(int j=0;j<n;j++)
	{
		f[i][j]=1;
	}
	row=m;
	col=n;
	return TRUE;
}
template<class Type>
BOOL CMatrix<Type>::Zeros(int m)
{ 	
	fMatrixAlloc(m,m);
	for(int i=0;i<m;i++)
	for(int j=0;j<m;j++)
	{
		f[i][j]=0;
	}
	row=m;
	col=m;
	return TRUE;
}
template<class Type>
BOOL CMatrix<Type>::Zeros(int m,int n)
{ 	
	fMatrixAlloc(m,n);
	if(f==NULL) return FALSE;
	for(int i=0;i<m;i++)
	for(int j=0;j<n;j++)
	{
		f[i][j]=0;
	}
	row=m;
	col=n;
	return TRUE;
}

template<class Type>
CMatrix<Type> CMatrix<Type>::operator*(CMatrix<Type>  &mat)
{
	int l,m,n;
	l=GetRow();
	m=GetCol();
	n=mat.GetCol();
	CMatrix tmp(l,n);
	tmp.Zeros(l,n);
	for(int k=0;k<l;k++)
	for(int j=0;j<n;j++)
	{
		for(int i=0;i<m;i++)
		tmp.f[k][j]+=f[k][i]*mat.f[i][j];
	}	
	return tmp;
}

template<class Type>
CMatrix<Type> CMatrix<Type>::operator/(Type x)
{
	int l,n;
	l=GetRow();
	n=GetCol();
	CMatrix<Type> tmp(l,n);
	for(int j=0;j<l;j++)
	for(int i=0;i<n;i++)
		tmp.f[j][i]=f[j][i]/x;
	return tmp;
}

template<class Type>
CMatrix<Type> CMatrix<Type>::operator+( CMatrix<Type>  &mat)
{
	int l1,m1;
	l1=GetRow();
	m1=GetCol();
	CMatrix<Type> tmp(l1,m1);
	for(int j=0;j<l1;j++)
	for(int i=0;i<m1;i++)
		tmp.f[j][i]=f[j][i]+mat.f[j][i];
	return tmp;
}

template<class Type>
CMatrix<Type> CMatrix<Type>::operator+(Type x)
{
	int l,m;
	l=GetRow();
	m=GetCol();
	CMatrix<Type> tmp(l,m);
	for(int j=0;j<l;j++)
	for(int i=0;i<m;i++)
		tmp.f[j][i]=f[j][i]+x;
	return tmp;
}

template<class Type>
CMatrix<Type> CMatrix<Type>::operator-(CMatrix<Type>  &mat)
{
	int l1,m1;
	l1=GetRow();
	m1=GetCol();
	CMatrix<Type> tmp(l1,m1);
	for(int j=0;j<l1;j++)
	for(int i=0;i<m1;i++)
		tmp.f[j][i]=f[j][i]-mat.f[j][i];
	return tmp;
}
template<class Type>
CMatrix<Type> CMatrix<Type>::operator-()
{
	for(int i=0;i<row;i++)
	for(int j=0;j<col;j++)
		f[i][j]=-f[i][j];
	return *this;
}	

template<class Type>
CMatrix<Type> CMatrix<Type>::operator-( Type x)
{
	int l,m;
	l=GetRow();
	m=GetCol();
	CMatrix<Type> temp(l,m);
	for(int j=0;j<l;j++)
	for(int i=0;i<m;i++)
		temp.f[j][i]=f[j][i]-x;
	return temp;
}

template<class Type>
istream &operator>>(istream &in,CMatrix<Type> &mat)
{    
	int i,j;
	for(i=0;i<mat.row;i++)
	for(j=0;j<mat.col;j++)
		 in>>mat.f[i][j];
	cout<<endl;
	return in;
};

template<class Type>
ostream  & operator<<(ostream& out,CMatrix<Type> &v1)
{
	if(v1.GetF()==NULL) 
	{
		out<<"This Matrix cannot be output!";
		return out<<endl;
	}
	out<<setiosflags(ios::right||ios::fixed);
	out<<setprecision(5); 
	out<<"["<<endl;
	int mr=v1.GetRow();
	int mc=v1.GetCol();
	for(int j=0;j<mr;j++)
	{
		for(int i=0;i<mc-1;i++)
			out<<setw(12)<<v1.f[j][i];
		out<<setw(12)<<v1.f[j][mc-1]<<";"<<endl;
	}
	return out<<"]"<<endl;
	
}
template<class Type>
CMatrix<Type>::CMatrix(CMatrix<Type> &copy)
{
	f=NULL;
	*this=copy;
}
	
template<class Type>
void CMatrix<Type>::operator=(const CMatrix<Type> &copy)
{
	
	if(this==&copy) return;
	if(f!=NULL) fMatrixFree();
	int m,n;
	m=copy.row;
	n=copy.col;
	fMatrixAlloc(m,n);
  	for(int i=0;i<m;i++)
	for(int j=0;j<n;j++)
		f[i][j]=copy.f[i][j];
	row=m;
	col=n;  
}

template<class Type>
void CMatrix<Type>::operator=(Type d)
{
	
	if(f==NULL) 
	{
		cout<<"The Matrix is not be allociated"<<endl;
		return;
	}
	for(int i=0;i<row;i++)
	for(int j=0;j<col;j++)
		f[i][j]=d;
}

template<class Type>
Type ** CMatrix<Type>::MatrixAlloc(int r,int c)
{
    Type *x,**y;
    int   n;
    x=(Type *)calloc(r*c,sizeof(Type));
    y=(Type **)calloc(r,sizeof(Type *));
	for(n=0;n<=r-1;++n)
	y[n]=&x[c*n];
	free(x);
	free(y);
	return(y);
}

template<class Type>
void CMatrix<Type>::fMatrixAlloc(int r,int c)
{
    f = new Type *[r];
	for(int i =0 ; i < r;i++)
	{
		f[i] = new Type[c];
	}
}

template<class Type>
void CMatrix<Type>::MatrixFree(Type **x)
{
   free(x[0]);
   free(x);
}

template<class Type>
void CMatrix<Type>::fMatrixFree()
{
	for(int i =0 ; i < row;i++)
	{
		delete []f[i];
	}
	delete []f;
}

template<class Type>
BOOL CMatrix<Type>::Malloc(int m,int n)
{
	fMatrixAlloc(m,n);
	if(f==NULL) return FALSE;
	row=m;
	col=n;
	return TRUE;
}

	
template<class Type>
CMatrix<Type> CMatrix<Type>::Turn()
{
	CMatrix<Type> tmp(col,row);
	for(int i=0;i<row;i++)
	for(int j=0;j<col;j++)
		tmp.f[j][i]=f[i][j];
	return tmp;
}

template<class Type>
BOOL CMatrix<Type>::Diag(Type *array,int m)
{
	fMatrixAlloc(m,m);
	for(int i=0;i<m;i++)
	for(int j=0;j<m;j++)
		if(i==j) f[i][j]=array[i];
	row=m;
	col=m;
	return TRUE;
}

template<class Type>
CMatrix<Type> CMatrix<Type>::Outer(CMatrix x,CMatrix y) //x.size=y.size
{
	int row = x.GetRow();
	int col = x.GetCol();
	CMatrix<Type> tmp(row,col);
	for(int i = 0;i<row; i++)
	{
		for(int j=0; j < col;j++)
		{
			tmp[i][j] = x[i][j]*y[i][j];
		}
	}
	return tmp;

}

template<class Type>
CMatrix<Type> CMatrix<Type>::Construct(CMatrix<Type> &mat,CString style)
{
	int i,j;
	CMatrix<Type> tmp;
	if(style=="LR"||style=="lr")
	{
		if(row!=mat.row)    return tmp;
		if(!tmp.Malloc(row,col+mat.col)) return tmp;
		for(i=0;i<tmp.row;i++)
		{
			for(j=0;j<tmp.col;j++)
			{
				if(j<col) tmp.f[i][j]=f[i][j];
				else      tmp.f[i][j]=mat.f[i][j-col];
			}
		}
	}
	return tmp;
}

template<class Type>	
CMatrix<Type> CMatrix<Type>::GetChild(int sr,int sc,int er,int ec)
{
	int i,j;
	CMatrix<Type> tmp(er-sr+1,ec-sc+1);
	for(i=0;i<tmp.row;i++)
	for(j=0;j<tmp.col;j++)
		tmp[i][j]=f[sr+i][sc+j];
	return tmp;
}
template<class Type>	
BOOL CMatrix<Type>::IsInverse()
{  
	CMatrix<Type> tmp(col,row);
	Type d,p;
    int *is,*js,i,j,k,v;
    if(row!=col)  
	{  
		//cout<<"\nrow!=column,this matrix can't be inversed";
		return false;    
	}
	for(i=0;i<tmp.row;i++)
	for(j=0;j<tmp.col;j++)
		tmp[i][j]=f[i][j];
    is=new int[row];
    js=new int[row];
    for(k=0;k<=row-1;k++)
	{
		d=0.0;
	    for(i=k;i<=row-1;i++)
	    for(j=k;j<=row-1;j++)
	    {  
			p=fabs(tmp.f[i][j]);
		    if(p>d) 
			{  
				d=p;
				is[k]=i;
				js[k]=j;  
			}
		}
	    if(d+1.0==1.0)     //singular
	    {  
			delete is,js;
		   // cerr<<"singular,this matrix can't be inversed\n";
		    return false;    
		}
	    if(is[k]!=k)
	    for(j=0;j<=row-1;j++)
		{   
			v=is[k];
		    p=tmp.f[k][j];
			tmp.f[k][j]=tmp.f[v][j];
			tmp.f[v][j]=p;
		}
	    if(js[k]!=k)
	    for(i=0;i<=row-1;i++)
		{  
			v=js[k];
		    p=tmp.f[i][k]; 
			tmp.f[i][k]=tmp.f[i][v];
			tmp.f[i][v]=p;
		}
	    tmp.f[k][k]=1.0/tmp.f[k][k];
	    for(j=0;j<=row-1;j++)
	    if(j!=k)
	       tmp.f[k][j]*=tmp.f[k][k];
	    for(i=0;i<=row-1;i++)
			if(i!=k)
	        for(j=0;j<=row-1;j++)
				if(j!=k)
		             tmp.f[i][j]-=tmp.f[i][k]*tmp.f[k][j];
		for(i=0;i<=row-1;i++)
	    if(i!=k)
			 tmp.f[i][k]=-tmp.f[i][k]*tmp.f[k][k];
	}
    // change row and column after inverse
	for(k=row-1;k>=0;k--)
	{  
		if(js[k]!=k)
		for(j=0;j<=row-1;j++)
		{  
			v=js[k];
		    p=tmp.f[k][j];
			tmp.f[k][j]=tmp.f[v][j];
			tmp.f[v][j]=p;
		}
	    if(is[k]!=k)
			 for(i=0;i<=row-1;i++)
			 {  
				 v=is[k];
				 p=tmp.f[i][k];  
				 tmp.f[i][k]=tmp.f[i][v];
				 tmp.f[i][v]=p;
		     }
      }
	delete is,js;
	return true;    
} 

template<class Type>	
CMatrix<double> CMatrix<Type>::Inverse()
{  
	CMatrix<Type> tmp(col,row);
	Type d,p;
    int *is,*js,i,j,k,v;
    if(row!=col)  
	{  
		cout<<"\nrow!=column,this matrix can't be inversed";
		return *this;    
	}
	for(i=0;i<tmp.row;i++)
	for(j=0;j<tmp.col;j++)
		tmp[i][j]=f[i][j];
    is=new int[row];
    js=new int[row];
    for(k=0;k<=row-1;k++)
	{
		d=0.0;
	    for(i=k;i<=row-1;i++)
	    for(j=k;j<=row-1;j++)
	    {  
			p=fabs(tmp.f[i][j]);
		    if(p>d) 
			{  
				d=p;
				is[k]=i;
				js[k]=j;  
			}
		}
	    if(d+1.0==1.0)     //singular
	    {  
			delete is,js;
		    cerr<<"\nerror*****not inv,be careful your matrix had been changed\n";
		    return *this;    
		}
	    if(is[k]!=k)
	    for(j=0;j<=row-1;j++)
		{   
			v=is[k];
		    p=tmp.f[k][j];
			tmp.f[k][j]=tmp.f[v][j];
			tmp.f[v][j]=p;
		}
	    if(js[k]!=k)
	    for(i=0;i<=row-1;i++)
		{  
			v=js[k];
		    p=tmp.f[i][k]; 
			tmp.f[i][k]=tmp.f[i][v];
			tmp.f[i][v]=p;
		}
	    tmp.f[k][k]=1.0/tmp.f[k][k];
	    for(j=0;j<=row-1;j++)
	    if(j!=k)
	       tmp.f[k][j]*=tmp.f[k][k];
	    for(i=0;i<=row-1;i++)
			if(i!=k)
	        for(j=0;j<=row-1;j++)
				if(j!=k)
		             tmp.f[i][j]-=tmp.f[i][k]*tmp.f[k][j];
		for(i=0;i<=row-1;i++)
	    if(i!=k)
			 tmp.f[i][k]=-tmp.f[i][k]*tmp.f[k][k];
	}
    // change row and column after inverse
	for(k=row-1;k>=0;k--)
	{  
		if(js[k]!=k)
		for(j=0;j<=row-1;j++)
		{  
			v=js[k];
		    p=tmp.f[k][j];
			tmp.f[k][j]=tmp.f[v][j];
			tmp.f[v][j]=p;
		}
	    if(is[k]!=k)
			 for(i=0;i<=row-1;i++)
			 {  
				 v=is[k];
				 p=tmp.f[i][k];  
				 tmp.f[i][k]=tmp.f[i][v];
				 tmp.f[i][v]=p;
		     }
      }
	delete is,js;
	return tmp;    
} 

template<class Type>	
CMatrix<Type> operator * (double num,CMatrix<Type> right)
{
    CMatrix<Type> tmp(right.row,right.col);
	int i,j;
	for(i=0;i<right.row;i++)
	for(j=0;j<right.col;j++)
		 tmp.f[i][j]=num*right.f[i][j];
	return tmp;
}

template<class Type>	
CMatrix<Type> operator * (CMatrix<Type> left,double num)
{    
	CMatrix<Type> tmp(left.row,left.col);
	int i,j;
	for(i=0;i<left.row;i++)
	for(j=0;j<left.col;j++)
		 tmp.f[i][j]=num*left.f[i][j];
	return tmp;
}
template<class Type>	
void CMatrix<Type>::Hessenberg()
{
	int i,j,k;
	double d,t;
	if(row==0)
	{
	    cerr<<"\na null matrix,can't use hessenberg transformation";
	    exit(1);
	 }
	if(row!=col)
	{
	    cerr<<"\nnot a square matrix,can't use hessenberg transformation";
	    exit(1);
	 }
	for(k=1;k<=row-2;k++)
	{
	    d=0.0;
	    for(j=k;j<=row-1;j++)
	    {
		 t=f[j][k-1];
		 if(fabs(t)>fabs(d))      { d=t;  i=j;  }
	     }
	    if(fabs(d)+1.0!=1.0)
	    {
		 if(i!=k)
		 {      for(j=k-1;j<=row-1;j++)
			{     t=f[i][j];  f[i][j]=f[k][j];  f[k][j]=t;
			 }
			for(j=0;j<=row-1;j++)
			{      t=f[j][i];  f[j][i]=f[j][k];  f[j][k]=t;
			 }
		  }
		 for(i=k+1;i<=row-1;i++)
		 {       t=f[i][k-1]/d;
			 f[i][k-1]=0.0;
			 for(j=k;j<=row-1;j++)
			 {      f[i][j]-=t*f[k][j];
			  }
			 for(j=0;j<=row-1;j++)
			 {      f[j][k]+=t*f[j][i];
			 }
		  }
	     }
	}
 }

template<class Type>	
BOOL CMatrix<Type>::Cholesky()
{
	int i,j,k;
	double d;
	if(row!=col)
	{	cerr<<"\nnot a squre matrix, can't use Cholesky disintegration";
		return FALSE;
	}
	if((f[0][0]+1.0==1.0)||(f[0][0]<0.0))
	{	cerr<<"\nnot a Zhengdin matrix, can't use Cholesky disintegration";
		return FALSE;
	}
	f[0][0]=sqrt(f[0][0]);
	d=f[0][0];
	for(i=1;i<=row-1;i++)
	{	f[i][0]/=f[0][0];
	}
	for(j=1;j<=row-1;j++)
	{	for(k=0;k<=j-1;k++)
		{   f[j][j]-=f[j][k]*f[j][k];
		}
        if((f[j][j]+1.0==1.0)||(f[j][j]<0.0))
		{   cerr<<"\nnot a zhendin matrix, can't use Cholesky disintegration";
		    return FALSE;
		}
		f[j][j]=sqrt(f[j][j]);
		d*=f[j][j];
	    for(i=j+1;i<=row-1;i++)
		{	for(k=0;k<=j-1;k++)
                f[i][j]-=f[i][k]*f[j][k];
		    f[i][j]/=f[j][j];
		}
	}
	for(i=0;i<=row-2;i++)
		for(j=i+1;j<=row-1;j++)
			f[i][j]=0.0;
    return TRUE;
}
template<class Type>	
CMatrix<Type> CMatrix<Type>::Jacobi(double eps,int jt)
{ 
		CMatrix<Type> v;
		v.Zeros(row,col);
		if(row!=col) return v;
		int i,j,p,q,l,n;
		n=row;
		for(i=0;i<=n-1;i++)    // check the matrix's symmetrition
			for(j=0;j<=n-1;j++)
				if((f[i][j]-f[j][i])>0.001) 
					return v;
		double fm,cn,sn,omega,x,y,d;
		l=1;
		for (i=0; i<=n-1; i++)
		{ 
			v[i][i]=1.0;
			for (j=0; j<=n-1; j++)
			if (i!=j) v[i][j]=0.0;
		}
    while (1==1)//????????????????????????????????????
      { fm=0.0;
        for (i=0; i<=n-1; i++)
        for (j=0; j<=n-1; j++)
          { d=fabs(f[i][j]);
            if ((i!=j)&&(d>fm))
              { fm=d; p=i; q=j;}
          }
        if (fm<eps)  return v;
        if (l>jt)  return v;
        l=l+1; 
        x=-f[p][q]; y=(f[q][q]-f[p][p])/2.0;
        omega=x/sqrt(x*x+y*y);
        if (y<0.0) omega=-omega;
        sn=1.0+sqrt(1.0-omega*omega);
        sn=omega/sqrt(2.0*sn);
        cn=sqrt(1.0-sn*sn);
        fm=f[p][p];
        f[p][p]=fm*cn*cn+f[q][q]*sn*sn+f[p][q]*omega;
        f[q][q]=fm*sn*sn+f[q][q]*cn*cn-f[p][q]*omega;
        f[p][q]=0.0; f[q][p]=0.0;
        for (j=0; j<=n-1; j++)
        if ((j!=p)&&(j!=q))
        {
            fm=f[p][j];
            f[p][j]=fm*cn+f[q][j]*sn;
            f[q][j]=-fm*sn+f[q][j]*cn;
          }
        for (i=0; i<=n-1; i++)
          if ((i!=p)&&(i!=q))
          {
              fm=f[i][p];
              f[i][p]=fm*cn+f[i][q]*sn;
              f[i][q]=-fm*sn+f[i][q]*cn;
            }
        for (i=0; i<=n-1; i++)
          {
            fm=v[i][p];
            v[i][p]=fm*cn+v[i][q]*sn;
            v[i][q]=-fm*sn+v[i][q]*cn;
          }
      }
    return v;
}


template<class Type>
int CMatrix<Type>::rank()
{ 
	int i,j,k,nn,is,js,l,u,v,m,n,ll;
	double q,d;
	Type *a;
	m = row;
	n = col;
	nn=row;
	if (row>=col) nn=col;
	a=(Type*)calloc(m*n,sizeof(Type));
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			a[i*n+j]=f[i][j];
		}
	k=0;
    for (l=0; l<=nn-1; l++)
      { q=0.0;
        for (i=l; i<=m-1; i++)
        for (j=l; j<=n-1; j++)
          { ll=i*n+j; d=fabs(a[ll]);
	    if (d>q) { q=d; is=i; js=j;}
          }
        if (q+1.0==1.0) return(k);
        k=k+1;
        if (is!=l)
          { for (j=l; j<=n-1; j++)
              { u=l*n+j; v=is*n+j;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        if (js!=l)
          { for (i=l; i<=m-1; i++)
              { u=i*n+js; v=i*n+l;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        ll=l*n+l;
        for (i=l+1; i<=n-1; i++)
          { d=a[i*n+l]/a[ll];
            for (j=l+1; j<=n-1; j++)
              { u=i*n+j;
                a[u]=a[u]-d*a[l*n+j];
              }
          }
      }
	free(a);
    return(k);
}

template<class Type>
CMatrix<double> CMatrix<Type>::QREigen(double eps,int jt)
{              // return eigenvalue in two row matrix of a real matrix
	int m,it,i,j,k,l,ii,jj,kk,ll;
    CMatrix<double> uv;
	uv.Zeros(row,2);
	if(row!=col) return uv;
	int n=row;
    double b,c,w,g,xy,p,q,r,x,s,e,ff,z,y;
    double *u,*v,*a;
    u=(double*)calloc(n,sizeof(double));
    v=(double*)calloc(n,sizeof(double));
    a=(double*)calloc(n*n,sizeof(double));
    this->Hessenberg();
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			a[i*n+j]=f[i][j];
	it=0; m=n;
    while (m!=0)
      { l=m-1;
        while ((l>0)&&(fabs(a[l*n+l-1])>eps*
	      (fabs(a[(l-1)*n+l-1])+fabs(a[l*n+l])))) l=l-1;
        ii=(m-1)*n+m-1; jj=(m-1)*n+m-2;
        kk=(m-2)*n+m-1; ll=(m-2)*n+m-2;
        if (l==m-1)
          { u[m-1]=a[(m-1)*n+m-1]; v[m-1]=0.0;
            m=m-1; it=0;
          }
        else if (l==m-2)
          { b=-(a[ii]+a[ll]);
            c=a[ii]*a[ll]-a[jj]*a[kk];
            w=b*b-4.0*c;
            y=sqrt(fabs(w));
            if (w>0.0)
              { xy=1.0;
                if (b<0.0) xy=-1.0;
                u[m-1]=(-b-xy*y)/2.0;
                u[m-2]=c/u[m-1];
                v[m-1]=0.0; v[m-2]=0.0;
              }
            else
              { u[m-1]=-b/2.0; u[m-2]=u[m-1];
                v[m-1]=y/2.0; v[m-2]=-v[m-1];
              }
            m=m-2; it=0;
          }
        else
          { if (it>=jt)
              { cout<<"fail\n";
				free(u);
				free(v);
				free(a);
                return uv;
              }
            it=it+1;
            for (j=l+2; j<=m-1; j++)
              a[j*n+j-2]=0.0;
            for (j=l+3; j<=m-1; j++)
              a[j*n+j-3]=0.0;
            for (k=l; k<=m-2; k++)
              { if (k!=l)
                  { p=a[k*n+k-1]; q=a[(k+1)*n+k-1];
                    r=0.0;
                    if (k!=m-2) r=a[(k+2)*n+k-1];
                  }
                else
                  { x=a[ii]+a[ll];
                    y=a[ll]*a[ii]-a[kk]*a[jj];
                    ii=l*n+l; jj=l*n+l+1;
                    kk=(l+1)*n+l; ll=(l+1)*n+l+1;
                    p=a[ii]*(a[ii]-x)+a[jj]*a[kk]+y;
                    q=a[kk]*(a[ii]+a[ll]-x);
                    r=a[kk]*a[(l+2)*n+l+1];
                  }
                if ((fabs(p)+fabs(q)+fabs(r))!=0.0)
                  { xy=1.0;
                    if (p<0.0) xy=-1.0;
                    s=xy*sqrt(p*p+q*q+r*r);
                    if (k!=l) a[k*n+k-1]=-s;
                    e=-q/s; ff=-r/s; x=-p/s;
                    y=-x-ff*r/(p+s);
                    g=e*r/(p+s);
                    z=-x-e*q/(p+s);
                    for (j=k; j<=m-1; j++)
                      { ii=k*n+j; jj=(k+1)*n+j;
                        p=x*a[ii]+e*a[jj];
                        q=e*a[ii]+y*a[jj];
                        r=ff*a[ii]+g*a[jj];
                        if (k!=m-2)
                          { kk=(k+2)*n+j;
                            p=p+ff*a[kk];
                            q=q+g*a[kk];
                            r=r+z*a[kk]; a[kk]=r;
                          }
                        a[jj]=q; a[ii]=p;
                      }
                    j=k+3;
                    if (j>=m-1) j=m-1;
                    for (i=l; i<=j; i++)
                      { ii=i*n+k; jj=i*n+k+1;
                        p=x*a[ii]+e*a[jj];
                        q=e*a[ii]+y*a[jj];
                        r=ff*a[ii]+g*a[jj];
                        if (k!=m-2)
                          { kk=i*n+k+2;
                            p=p+ff*a[kk];
                            q=q+g*a[kk];
                            r=r+z*a[kk]; a[kk]=r;
                          }
                        a[jj]=q; a[ii]=p;
                      }
                  }
              }
          }
      }

    for(i=0;i<n;i++)
		uv[i][0]=u[i];
	for(i=0;i<n;i++)
		uv[i][1]=v[i];
	free(u);
	free(v);
	free(a);
	return uv;
  }
/*
template<class Type>
int CMatrix<Type>::Svd(Type *u, Type *v, double eps)
{
	int m,n,ka,i,j,flag;
	m = row;
	n = col;
	ka = (m>n?(m+1):(n+1));
	Type *a;
	a=(Type*)calloc(m*n,sizeof(Type));

	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			a[i*n+j]=f[i][j];
		}
	flag = dluav(a,m,n,u,v,eps,ka);
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			f[i][j] = a[i*n+j];
		}
	free(a);
	return flag;
}
*/
#endif

 
