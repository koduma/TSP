#pragma warning(disable:4710)
#pragma warning(disable:4711)
#pragma warning(disable:4820)
#pragma GCC target ("sse4.2")
#include <vector>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <climits>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <random>
#include <queue>
#include <deque>
#include <list>
#include <map>
#include <array>
#include <chrono>
#include <fstream>
#include <functional>
#include <unordered_map>
//#include <immintrin.h>
//#include "hash_map.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#define DIR 2392
#define CITY 2392
#define TH 2
    
typedef unsigned long long ull;

#define INF 10000000

ull zoblish_field[CITY][CITY];

struct node {//どういう手かの構造体
	int route[CITY+1];//スワイプ移動座標
	int pos;
	int t;
	int score;//評価値
	int score2;
	ull hash;//盤面のハッシュ値
	bool visited[CITY]={0};
	node() {//初期化
		this->score = 0;
	}
	bool operator < (const node& n)const {//スコアが高い方が優先される
		return score < n.score;
	}
}sim[TH];

ull xor128() {//xorshift整数乱数
	static unsigned long long rx = 123456789, ry = 362436069, rz = 521288629, rw = 88675123;
	ull rt = (rx ^ (rx << 11));
	rx = ry; ry = rz; rz = rw;
	return (rw = (rw ^ (rw >> 19)) ^ (rt ^ (rt >> 8)));
}
int rnd(int mini, int maxi) {
	static mt19937 mt((int)time(0));
	uniform_int_distribution<int> dice(mini, maxi);
	return dice(mt);
}
double d_rnd() {
	static mt19937 mt((int)time(0));
	uniform_real_distribution<double> dice(0, 1.0);
	return dice(mt);
}
node solve(node travel,int kyori[CITY][CITY],double checksum[CITY][CITY]) {
    
    int ttt=(travel.t)+1;
    
    int MLEN=INF;
    int P_MLEN=INF;
    
    node no1;
    
    unordered_map<ull,bool>mp;
       
    int discover=0;
    
    int counter=0;
    
    vector<node>elite;
    
    int iter=0;
    
    while(1){
        
    if(travel.score<=440000&&(int)elite.size()<=100){elite.push_back(travel);}    
        
    int choice=-1;    
	    
    if((int)elite.size()>0){        
    int r=rnd(0,(int)elite.size()-1);    
    travel=elite[r];
    choice=r;    
    }    
        
    node old=travel;
    node news=travel;
     
    for(int u=0;u<1;u++){    
        
    int alpha;
    if(ttt>CITY-3){alpha=CITY-3;}
    else{alpha=rnd(ttt,CITY-3);}
    int beta=rnd(alpha+1,CITY-1);
        
    int a=(int)news.route[alpha-1];
    int b=(int)news.route[alpha];
    int c=(int)news.route[beta];
    int d=(int)news.route[beta+1];    
        
    int prev=kyori[a][b]+kyori[c][d];
    int neww=kyori[a][c]+kyori[b][d];
        
    news.score+=neww-prev;
    vector<int>rrr;
    for(int i=beta;i>=alpha;i--){
    rrr.push_back(news.route[i]);
    }
    for(int i=0;i<(int)rrr.size();i++){
    news.route[alpha+i]=rrr[i];
    }    
        
    }    
        
    double nscore=(double)(old.score/news.score);//neww<prev->good
        
    int drnd=d_rnd();
        
    if(news.score<old.score || abs(old.score-news.score)<=100){ 
    //if(news.score<old.score){    
    discover++;
    travel=news;
    if(choice>=0){
    elite[choice]=news;
    }    
    }   
    if(MLEN>travel.score){
    MLEN=travel.score;
    no1=travel;
    }
    iter++; 
        
    if(iter%10000==0){
        P_MLEN=MLEN;
    }
    //if(iter%10000==0){cout<<iter/10000<<"%,score="<<(MLEN)<<endl;}    
    if(iter%1000000==0){break;}
    }
return no1;
}
void check_travel(node ans,double checksum[CITY][CITY]){
	
	int i;
	cout<<"ans.score="<<ans.score<<endl;
	int check[CITY]={0};

	double check_score=0;

	int pos=0;

	for(i=0;i<=CITY;i++){
		cout<<ans.route[i]+1<<",";
		check[ans.route[i]]++;
		check_score+=checksum[pos][(int)ans.route[i]];
		pos=(int)ans.route[i];
	}
	printf("\n");

	int ok=0;

	if(check[0]==2){ok++;}

	for(i=1;i<CITY;i++){
		if(check[i]==1){ok++;}
	}

	printf("ok=%d\n",ok);
	printf("check_score=%lf\n",check_score);

}

int XXX[CITY][CITY];
double YYY[CITY][CITY];

void sub(int kyori[CITY][CITY],double checksum[CITY][CITY]){
	
	FILE *fp=fopen("pr2392.txt","r");
	
	int i,j,k;
	double u[CITY][2];
	for(i=0;i<CITY;i++){
	int n;
	fscanf(fp,"%d %lf %lf",&n,&u[i][0],&u[i][1]);
	}
    
	
	for(i=0;i<CITY;i++){
	for(j=0;j<CITY;j++){
	kyori[i][j]=(int)ceil(sqrt((u[i][0]-u[j][0])*(u[i][0]-u[j][0])+(u[i][1]-u[j][1])*(u[i][1]-u[j][1])));
	checksum[i][j]=sqrt((u[i][0]-u[j][0])*(u[i][0]-u[j][0])+(u[i][1]-u[j][1])*(u[i][1]-u[j][1]));
	}
	}
	
	for(i=0;i<CITY;i++){
	for(j=0;j<CITY;j++){
	zoblish_field[i][j]=xor128();
	}
	}
    
	
	node travel,ans;

	travel.route[0]=0;
	travel.pos=0;
	travel.t=0;
	travel.score=0;
	for(i=0;i<CITY;i++){
		travel.visited[i]=false;
	}
	travel.visited[0]=true;
	
	int history[CITY+1];
        history[0]=1;
        history[CITY]=1;
        for(i=1;i<=CITY-1;i++){history[i]=i+1;}
	
	vector<node>dque;
	for(int x=0;x<TH;x++){
    for (i = 0; i < 100; i++) {
        int r1 = rnd(1,CITY-1);
        int r2 = rnd(1,CITY-1);
        int tmp = history[r2];
        history[r2] = history[r1];
        history[r1] = tmp;
    }    
        travel.score=0;
        for(i=0;i<CITY;i++){
            travel.route[i]=history[i]-1;
            travel.score+=kyori[history[i]-1][history[i+1]-1];
        }
        travel.route[CITY]=0;    
	dque.push_back(travel);
	}
	
	for (int ii = 0; ii < 1000; ii++) {
	int ks=(int)dque.size();
#pragma omp parallel for
	for(int x=0;x<ks;x++){
	node cand=dque[x];
	sim[x]=solve(cand,kyori,checksum);
	}
	int d=INF;
	int index=0;
	dque.clear();
	for (int x = 0; x < TH;x++) {
	if(d>sim[x].score){d=sim[x].score;index=x;}
	dque.push_back(sim[x]);
	}
	check_travel(sim[index],checksum);    
        }
}

int main(){
    
  int i,j,k;
    
  sub(XXX,YYY);

  cin>>i;
  cin>>j;
	cin>>k;
  
  return 0;
}
