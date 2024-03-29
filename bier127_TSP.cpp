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

#define BEAM_WIDTH 1
#define BEAM_WIDTH2 1
#define DIR 127
#define CITY 127
#define NODE_SIZE BEAM_WIDTH*DIR
#define AVG 3

typedef unsigned long long ull;

#define INF 10000000

int kyori[CITY][CITY];
double checksum[CITY][CITY];

int call=0;

int MLEN=INF;

ull zoblish_field[CITY][CITY];

struct node {//どういう手かの構造体
	char route[CITY+1];//スワイプ移動座標
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
}fff[NODE_SIZE],ff[DIR*BEAM_WIDTH2],no1;

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
node solve(node travel,int count) {
    
    

	vector<node>dque;

	dque.push_back(travel);
	//emilib::HashMap<ull, bool> checkNodeList;
	
	int ttt=(travel.t)+1;
    
    /*

	//2手目以降をビームサーチで探索
	for (int i = (travel.t)+1; i <= CITY; i++) {
		int ks = (int)dque.size();
#pragma omp parallel for
		for (int k = 0; k < ks; k++) {
			node temp = dque[k];//que.front(); que.pop();
			if(i==CITY){
				node cand = temp;
				int sc=cand.score;
				cand.score=sc+kyori[cand.pos][0];
				cand.route[i]=0;
				cand.pos=0;
				cand.t=i;
				cand.visited[0]=true;
				//cand.hash=calc_hash(cand);
				fff[k] = cand;
			}
			else{
			for (int j = 0; j < DIR; j++) {//上下左右の4方向が発生
				node cand = temp;
				if (!cand.visited[j]) {
					int sc=cand.score;
					cand.score=sc+kyori[cand.pos][j];
					cand.route[i]=j;
					cand.pos=j;
					cand.t=i;
					cand.visited[j]=true;
					//cand.hash=calc_hash(cand);
					fff[(CITY * k) + j] = cand;
						//part4 += omp_get_wtime() - st;
					}
					else {
						cand.score = -1;
						fff[(CITY * k) + j] = cand;
					}
			}
		}
		}
		//printf("depth=%d/%d\n",i,CITY);
		dque.clear();
		vector<pair<int,int> >vec;
		int ks2 = 0;
		if(i==CITY){
			for (int j = 0; j < ks; j++) {
				if (fff[j].score != -1) {
				vec.push_back(make_pair(fff[j].score,j));
				ks2++;
				}
			}
		}
		else{
		for (int j = 0; j < CITY * ks; j++) {
			if (fff[j].score != -1) {
			vec.push_back(make_pair(fff[j].score,j));
			ks2++;
			}
		}
	}

	sort(vec.begin(),vec.end());
	int push_node=0;
		for (int j = 0; push_node < BEAM_WIDTH && j<ks2 ;j++) {
			int v=vec[j].second;
			node temp = fff[v];
			if (i==CITY) {//コンボ数が増えたらその手を記憶する
				call++;
				if(call%10==0){
				//printf("call=%d/%d\n",call,CITY+(CITY*(CITY-1)/2)*BEAM_WIDTH2*AVG);
				}
				travel=temp;
				goto SA;
			}
			//if(temp.score>=5337){break;}
			//if(!checkNodeList[temp.hash]){
				//checkNodeList[temp.hash]=true;
				dque.push_back(temp);
				push_node++;
				//}
			}
	}
        
    }
    
    SA:
    */
    
    unordered_map<ull,bool>mp;
    
    mp.clear();
    
    if(count==0){
       
    int history[CITY+1]={1,16,2,51,50,115,13,120,10,100,64,58,91,61,62,59,60,116,90,3,11,9,24,23,4,22,19,72,8,67,73,74,68,71,70,69,75,76,78,117,84,81,126,82,83,101,102,63,119,96,109,88,87,86,85,110,104,125,89,92,99,65,113,66,55,124,52,5,56,121,57,54,44,103,45,48,118,53,47,49,46,94,112,111,107,127,93,95,123,97,98,32,29,28,122,33,25,26,38,39,42,34,43,40,35,37,36,41,14,12,30,27,31,80,79,77,18,21,17,20,108,15,106,6,114,105,7,1};
        travel.score=0;
        for(int i=0;i<CITY;i++){
            travel.route[i]=history[i]-1;
            travel.score+=kyori[history[i]-1][history[i+1]-1];
        }
        travel.route[CITY]=0;
    }    
    
    double start = omp_get_wtime();
    
    int discover=0;
    
    int counter=0;
    
    vector<node>elite;
    
    int iter=0;
    
    while(1){
    double elapsed = omp_get_wtime()-start;
    //if(elapsed>=10.0&&(int)elite.size()==0){break;}
    if(elapsed>=300.0){break;}
    
    if(5.0*(double)counter<=elapsed){cout<<5.0*(double)counter<<"s/"<<"300s"<<endl;counter++;cout<<"score="<<no1.score<<endl;}
        
    ull ha=0ll;
    
    for(int z=0;z<CITY;z++){
    ha^=zoblish_field[z][(int)travel.route[z]];
    }
        
    if(travel.score<121000){
        if(!mp[ha]){
        elite.push_back(travel);
        if((int)elite.size()%1000==0){
        cout<<"elite_size="<<(int)elite.size()<<endl;
        }
        mp[ha]=true;    
        }
    }
	    
    int choice=-1;
	    
    if((int)elite.size()>0){        
    int r=rnd(0,(int)elite.size()-1);    
    travel=elite[r];
    choice=r;
    //travel=no1;
    }
        
    node old=travel;
    node news=travel;
     
    for(int u=0;u<2;u++){    
        
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
    vector<char>rrr;
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
    discover++;
    travel=news;
    if(choice>=0){
    elite[choice]=news;
    }
    }   
    if(MLEN>travel.score){
    printf("score=%d\n",travel.score);
    MLEN=travel.score;
    no1=travel;
    }
    iter++;
    if(iter%10000000==0){cout<<"iter="<<iter<<endl;}
    }
	
return no1;
}
void check_travel(node ans){
	
	int i;
	cout<<"ans.score="<<ans.score<<endl;
	int check[CITY]={0};

	double check_score=0;

	int pos=0;

	for(i=0;i<=CITY;i++){
		cout<<(int)ans.route[i]+1<<",";
		check[(int)ans.route[i]]++;
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

int main(){  
	
	FILE *fp=fopen("bier127.txt","r");
	
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
	
	for(i=0;i<100;i++){
	ans=solve(travel,i);
	check_travel(no1);
	travel=ans;
	}
	
	cin>>i;
	cin>>j;
	cin>>k;

  return 0;
}
