#include <iostream>
#include <fstream> 
#include <string> 
#include <cstdlib> 
#include <sstream> 
#include <cmath> 
#include <ctime> 
#include <algorithm>

#define ROW 200 
#define COL 2
#define filename "Bubble.csv" 
using namespace std; 

double tmp;
const double PI  =3.141592653589793238463;
double Y[ROW][COL] = {0};
const int pallnum=42;
int pall[10000][pallnum];
int iter=0;
int itermax=10;
double tc,w,b,phi;

void read();
double fitLPPL(int i);
double LPPL(double A,double B,double C,double tc,double beta,double w,double phi);

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
 				Y[i][j] = atof(str_buf.c_str()); 
			} 
		} 
 		fin.close(); 
 		fin.clear(); 
	} 
}

double fitLPPL(int i){
	tmp=0;
	for(int j=0;j<8;j++){
		
		tmp+=pall[i][j]*(pow(2,j));
	}
	tc=200+tmp;
	
	tmp=0;
	for(int j=8;j<22;j++){
		tmp+=pall[i][j]*(pow(2,j-8));
	}
	w=0.01+0.01*tmp;
	
	tmp=0;
	for(int j=22;j<32;j++){
		tmp+=pall[i][j]*pow(2,j-22);
	}	
	b=(tmp+1)/1025;
	
	tmp=0;
	for(int j=32;j<42;j++){
		tmp+=pall[i][j]*pow(2,j-32);
	}
	phi=(tmp+1)/1025*PI*2;
	
	double M[199][3]={0};
	for(int i=0;i<199;i++){
		M[i][0]=1;
	}
	
	for(int i=0;i<199;i++){
		M[i][1]=pow(200-(i+1),b);
	}
	
	for(int i=0;i<199;i++){
		M[i][2]=pow((200-(i+1)),b)*cos(w*log(200-(i+1))+phi);
	}
	
	double N[199]={0};
	for(int i=0;i<199;i++){
		N[i]=log(Y[i][1]);
	}
	
	//Regress(N,M)
	double A,B,C,out;
	A=4.2746;
	B=-0.002367;
	C=0.0014;
	C=C/B;
	out=LPPL(A,B,C,tc,b,w,phi);
	
	double E;
	for(int i=0;i<199;i++){
		E=abs(out-N[i]);//out is array??mean?? 
	}
	return E;
}

double LPPL(double A,double B,double C,double tc,double beta,double w,double phi){
	for(int i=1;i<tc-1;i++){
		tc=200;
		return A+B*(pow(tc-i,beta))*(1+C*cos(w*log(tc-i)+phi));
	}
}

template<typename T> 
void bubble_sort(T arr[], int len) {
	int i, j;
	for (i = 0; i < len - 1; i++)
		for (j = 0; j < len - 1 - i; j++)
			if (arr[j][1] > arr[j + 1][1]){
				swap(arr[j][0], arr[j + 1][0]);
				swap(arr[j][1], arr[j + 1][1]);
			}			
}

int main(){
	read();
	
	srand(time(0));
	for(int i=0;i<10000;i++){
		for(int j=0;j<42;j++){
			pall[i][j]=rand()%2;
		}
	}
	
	double err[10000]={0};
	while(iter<10){
		for(int i=0;i<10000;i++){
			err[i]=fitLPPL(i);
		}
		
		double AA[10000][2]={0};
		for(int i=0;i<10000;i++){
			AA[i][0]=i;
			AA[i][1]=err[i];
		}
		
		int len = (int) sizeof(AA) / sizeof(*AA);
		bubble_sort(AA, len);
		
		cout<<iter<<" "<<fixed<<AA[0][1]<<endl;
		
		double psurvive[100][pallnum]={0};
		for(int i=0;i<100;i++){
			int tmp=AA[i][0];
			for(int j=0;j<pallnum;j++){
				psurvive[i][j]=pall[tmp][j];
			}	
		}
		
		for(int i=0;i<100;i++){
			for(int j=0;j<pallnum;j++){
				pall[i][j]=psurvive[i][j];
			}	
		}
		
		for(int i=100;i<10000;i++){
			int father=rand() % 100;//0~<100
			int mother=rand() % 100;
			int genepoint=rand() % pallnum-1;
			
			for(int j=0;j<genepoint;j++){
				pall[i][j]=pall[father][j];
			}
			
			for(int j=genepoint;j<pallnum;j++){
				pall[i][j]=pall[mother][j];
			}
		}
		
		for(int i=0;i<10000;i++){
			int son=(rand() %9900) +100;
			int genepoint=rand() % pallnum-1;
			pall[son][genepoint]=!pall[son][genepoint];
		}
		iter++;
	}
	
	for(int i=0;i<10000;i++){
		err[i]=fitLPPL(i);
	}
	
	//[tmp mini]=min(err); ???
	//code=pall(mini,:); ???
	int code[pallnum]={0};
	float min=100;
	int minnum;
	for(int i=0;i<10000;i++){
		if(err[i]<min){
			min=err[i];
			minnum=i;
		}
	}
	
	for(int i=0;i<pallnum;i++){
		code[i]=pall[minnum][i];
	}
	
	//resemble function "filLPPL"
	tmp=0;
	for(int j=0;j<8;j++){
		tmp+=code[j]*(pow(2,j));
	}
	tc=200+tmp;
	
	tmp=0;
	for(int j=8;j<22;j++){
		tmp+=code[j]*(pow(2,j-8));
	}
	w=0.01+0.01*tmp;
	
	tmp=0;
	for(int j=22;j<32;j++){
		tmp+=code[j]*pow(2,j-22);
	}	
	b=(tmp+1)/1025;
	
	tmp=0;
	for(int j=32;j<42;j++){
		tmp+=code[j]*pow(2,j-32);
	}
	phi=(tmp+1)/1025*PI*2;
	
	const int arrSize=tc-1;
	double M[arrSize][3]={0};
	for(int i=0;i<tc-1;i++){
		M[i][0]=1;
	}
	
	for(int i=0;i<tc-1;i++){
		M[i][1]=pow(tc-1-(i+1),b);
	}
	
	for(int i=0;i<tc-1;i++){
		M[i][2]=pow((tc-1-(i+1)),b)*cos(w*log(tc-1-(i+1))+phi);
	}
	
	double N[arrSize]={0};
	for(int i=0;i<tc-1;i++){
		N[i]=log(Y[i][1]);
	}
	
	//Regress(N,M)
	double A,B,C,out;
	A=4.2746;
	B=-0.002367;
	C=0.0014;
	C=C/B;
	out=LPPL(A,B,C,tc,b,w,phi);
	
	double E;
	for(int i=0;i<199;i++){
		E=abs(out-N[i]);//out is array??mean?? 
	}
	
	return 0;
}
