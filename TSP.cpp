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
#include <immintrin.h>
#include "hash_map.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#define BEAM_WIDTH 1000000
#define BEAM_WIDTH2 2
#define DIR 100
#define CITY 100
#define NODE_SIZE BEAM_WIDTH*DIR

typedef unsigned long long ull;

int kyori[CITY][CITY];

ull zoblish_field[CITY+1][CITY];

int call=0;

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
}fff[NODE_SIZE],ff[DIR*BEAM_WIDTH2];

ull calc_hash(node travel){

ull hash=0ll;

for (int i = 0; i <=travel.t; i++) {
char num = travel.route[i];
hash ^= zoblish_field[i][(int)num];
}

return hash;

}

ull xor128() {//xorshift整数乱数
	static unsigned long long rx = 123456789, ry = 362436069, rz = 521288629, rw = 88675123;
	ull rt = (rx ^ (rx << 11));
	rx = ry; ry = rz; rz = rw;
	return (rw = (rw ^ (rw >> 19)) ^ (rt ^ (rt >> 8)));
}

node BEAM_SEARCH(node travel) {

	vector<node>dque;

	dque.push_back(travel);
	emilib::HashMap<ull, bool> checkNodeList;

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
				cand.hash=calc_hash(cand);
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
					cand.hash=calc_hash(cand);
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
		printf("depth=%d/%d\n",i,CITY);
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
				printf("call=%d\n",call);
				return temp;
			}
			if(!checkNodeList[temp.hash]){
				checkNodeList[temp.hash]=true;
				dque.push_back(temp);
				push_node++;
				}
			}
	}

	return travel;
}
node BEAM_SEARCH2(node travel) {

	/*

	char route[CITY+1];//スワイプ移動座標
  int pos;
	int t;
	int score;//評価値
	ull hash;//盤面のハッシュ値
  bool visited[CITY]={0};

	*/


	vector<node>dque;

	dque.push_back(travel);
	emilib::HashMap<ull, bool> checkNodeList;

	//2手目以降をビームサーチで探索
	for (int i = (travel.t)+1; i <= CITY; i++) {
		int ks = (int)dque.size();
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
				cand.hash=calc_hash(cand);
				cand.score2=cand.score;
				ff[k] = cand;
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
					cand.hash=calc_hash(cand);
					node bbb=BEAM_SEARCH(cand);
					cand.score2=bbb.score;
					ff[(CITY * k) + j] = cand;
						//part4 += omp_get_wtime() - st;
					}
					else {
						cand.score2 = -1;
						cand.score = -1;
						ff[(CITY * k) + j] = cand;
					}
			}
		}
		}
		printf("depth2=%d/%d\n",i,CITY);
		dque.clear();
		vector<pair<int,int> >vec;
		int ks2 = 0;
		if(i==CITY){
			for (int j = 0; j < ks; j++) {
				if (ff[j].score != -1) {
				vec.push_back(make_pair(ff[j].score2,j));
				ks2++;
				}
			}
		}
		else{
		for (int j = 0; j < CITY * ks; j++) {
			if (ff[j].score != -1) {
			vec.push_back(make_pair(ff[j].score2,j));
			ks2++;
			}
		}
	}
	sort(vec.begin(),vec.end());
	int push_node=0;
		for (int j = 0; push_node < BEAM_WIDTH2 && j<ks2 ;j++) {
			int v=vec[j].second;
			node temp = ff[v];
			if (i==CITY) {//コンボ数が増えたらその手を記憶する
				return temp;
			}
			if(!checkNodeList[temp.hash]){
				checkNodeList[temp.hash]=true;
				dque.push_back(temp);
				push_node++;
				}
			}
	}

	return travel;
}

int main(){

  FILE *fp=fopen("kyori.txt","r");

  int i,j,k;

  for(i=0;i<CITY;i++){
    for(j=0;j<CITY;j++){
      fscanf(fp,"%d",&kyori[i][j]);
    }
  }

	for(i=0;i<CITY+1;i++){
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
	travel.hash=calc_hash(travel);

	ans=BEAM_SEARCH(travel);

	cout<<"ans.score="<<ans.score<<endl;

	int check[100]={0};

	int check_score=0;

	int pos=0;

	for(i=0;i<=CITY;i++){
		cout<<(int)ans.route[i]<<",";
		check[(int)ans.route[i]]++;
		check_score+=kyori[pos][(int)ans.route[i]];
		pos=(int)ans.route[i];
	}
	printf("\n");

	int ok=0;

	if(check[0]==2){ok++;}

	for(i=1;i<CITY;i++){
		if(check[i]==1){ok++;}
	}

	printf("ok=%d\n",ok);
	printf("check_score=%d\n",check_score);

	/*

	char route[CITY+1];//スワイプ移動座標
	int pos;
	int t;
	int score;//評価値
	ull hash;//盤面のハッシュ値
	bool visited[CITY]={0};

	*/


  cin>>i;
  cin>>j;
  cin>>k;

  return 0;
}
