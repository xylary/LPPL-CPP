#include <iostream>
#include <fstream> /* ifstream class */
#include <string> /* string class */
#include <cstdlib> /* atoi */
#include <sstream> /* istringstream class */
#include <cmath> 
#include <ctime> 
using namespace std; 
#define ROW 200 
#define COL 2
#define filename "Bubble.csv" 

const double PI  =3.141592653589793238463;
long double arr[ROW][COL] = {0}; 
int M[200-1][3]={0};
long double N[200]={0};
long double A=0;
long double B=0;
long double C=0;
long double tc[10000]={0};
long double w[10000]={0};
long double b[10000]={0};
long double phi[10000]={0};
long double out[10000]={0};
long double err[10000]={0};
long double E=0;
long double Aa[10000][2]={0};
long double Bb[10000][2]={0};
int test[42]= {1,0,1,0,0,0,1,1,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,1,0,1,0,0,0,0,0,0};
//rdsf 
int iter=0;
int itermax=10;

const int pallnum=42;
int pall[10000][pallnum];
int psurvive[100][42];

void read();
int fitLPPL();
void LPPL(long double A,long double B,long double C,long double tc,long double beta,long double w,long double phi);

void read(){
	ifstream fin(filename); 
	if(!fin) 
		cout << "open fail." << endl; 
	else{ 
		string str_buf; 
		int i=0, j=0; 
		
		for(i=0; i!=ROW; ++i){ 
			getline(fin, str_buf);
 	
			istringstream stream(str_buf); 
 
 			for(j=0; j!=COL; ++j) { 
 
	 			getline(stream, str_buf, ',');  
 				arr[i][j] = atof(str_buf.c_str()); 
			} 
		} 
 		fin.close(); 
 		fin.clear(); 
 		
		for(i=0; i!=ROW; ++i){ 
 			for(j=0; j!=COL; ++j) 
 				cout <<fixed<< arr[i][j] << ' '; 
 				cout << endl; 
		} 
	} 
}

int fitLPPL(){
	
	for(int i=0;i<199;i++){
		M[i][0]=1;
	}
	
	for(int i=0;i<199;i++){
		M[i][1]=pow(200-(i+1),b[i]);
	}
	
	for(int i=0;i<199;i++){
		M[i][2]=pow((200-(i+1)),b[i])*cos(w[i]*log(200-(i+1))+phi[i]);
	}
	
	for(int i=0;i<199;i++){
		N[i]=log(arr[i][1]);
	}
	
	
	//Regress(N,M)
	A=4.3;
	B=-2.36;
	C=0.0014;
	C=C/B;
	
	for(int i=0;i<199;i++){
		LPPL(A,B,C,tc[i],b[i],w[i],phi[i]);
	}//receive "out"
	for(int i=0;i<199;i++){
		out[i]=abs(out[i]-N[i]);
	}
	int tmp;
	for(int i=0;i<199;i++){
		tmp+=out[i];
	}
	E=tmp/199;
	return E;
}

void LPPL(long double A,long double B,long double C,long double tc,long double beta,long double w,long double phi){
	for(int i=1;i<tc-1;i++){
		out[i] = A+B*(pow(tc-i,beta))*(1+C*cos(w*log(tc-i)+phi));
	}
}

int main(){ 
	/*
	 1
     0
     1
     0
     0
     0
     1
     1
     0
     1
     1
     0
     1
     0
     0
     1
     0
     0
     0
     0
     0
     0
     0
     0
     0
     1
     0
     0
     0
     1
     1
     1
     0
     1
     0
     1
     0
     0
     0
     0
     0
     0

   397

    1.5100

    0.8829

    0.0674
	*/
	int i=0;
	float tmp=0;
	for(int j=0;j<42;j++){
		
		cout<<test[j];
	}
	cout<<"**************"<<endl;
	for(int j=0;j<8;j++){
		
		tmp+=test[j]*(pow(2,j));
	}
	tc[i]=200+tmp;
	cout<<tc[i]<<endl;
	//OK
	tmp=0;
	for(int j=8;j<22;j++){
		tmp+=test[j]*(pow(2,j-8));
	}
	w[i]=0.01+0.01*tmp;
	cout<<w[0]<<endl;
	//OK
	tmp=0;
	for(int j=22;j<32;j++){
		tmp+=test[j]*pow(2,j-22);
	}	
	b[i]=(tmp+1)/1025;
	cout<<b[0]<<endl;
	//OK
	tmp=0;
	for(int j=32;j<42;j++){
		tmp+=test[j]*pow(2,j-32);
	}
	phi[i]=(tmp+1)/1025*PI*2;
	cout<<phi[i]<<endl;
	//OK
	
	//Aa[i][1]=fitLPPL();	
	/*
	srand(time(0));
	read();
	
	for(int i=0;i<10000;i++){
		for(int j=0;j<42;j++){
			pall[i][j]=rand()%2;
		}
	}
	*/
	/*
	while(iter<itermax){
		for(int i=0;i<100;i++){//i<10000 
			int tmp=0;
			for(int j=0;j<8;j++){
				
				tmp+=pall[i][j]*(pow(2,j));
			}
			tc[i]=200+tmp;
			
			for(int j=8;j<22;j++){
				int tmp;
				tmp+=pall[i][j]*(pow(2,j-8));
			}
			w[i]=0.01+0.01*tmp;
	
			for(int j=22;j<32;j++){
				int tmp;
				tmp+=pall[i][j]*pow(2,j-22);
			}	
			b[i]=(tmp+1)/1025;
			
			for(int j=32;j<42;j++){
				int tmp;
				tmp+=pall[i][j]*pow(2,j-32);
			}
			phi[i]=(tmp+1)/1025*PI*2;
			
			Aa[i][1]=fitLPPL();	
	 }
		
		for(int i=0;i<10000;i++){
			Aa[i][0]=i;
			err[i];
		}
		*/
		/*
		long double k;
		for(int i=0;i<10000-1;i++){ //******
			if(Aa[i][1] > Aa[i+1][1]){ 
				k=Aa[i][1]; 
				Aa[i][1]=Aa[i-1][1]; 
				Aa[i-1][1]=k; 
			}
		}
				
		for(int i=0;i<10000;i++){
			//cout<<Aa[0][i]<<endl;
		}
		
		cout<<iter<<" "<<Aa[0][1]<<endl;
		
		for(int i=0;i<100;i++){
			for(int j=0;j<42;j++){
				psurvive[i][j]=pall[i][j];//???? 
			}
		}
		
		for(int i=0;i<100;i++){
			for(int j=0;j<42;j++){
				pall[i][j]=psurvive[i][j];
			}
		}
		
		for(int i=100;i<10000;i++){
			int father=rand() % 100;
			int mother=rand() % 100;
			int genepoint=rand() % 41;
			
			for(int j=0;j<genepoint;j++){
				pall[i][j]=pall[father][j];
			}
			
			for(int j=genepoint;j<42;j++){
				pall[i][j]=pall[mother][j];
			}
			
		}
		
		for(int i=1;i<10000;i++){
			int son=(rand() %9900) +100;
			int genepoint=rand() %41;
			pall[son][genepoint]=!pall[son][genepoint];
		}
		
		iter++;
	}
	*/
 	return 0;
} 
