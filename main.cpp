// main.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<random>
//#include<stdlib.h>
#include<time.h>
#include<unordered_map>
#include<algorithm>
#include"coor.h"
using namespace std;
//基本仿真参数设置
const int W = 20000; //带宽
const int H = 100;//UAV高度
const double v_max = 26.0; //UAV最大速度
const int alpha = 2; //路径损耗指数
const int beta = pow(10, 8);
const int node_num = 10; //待采集节点数目
const double err = 0.0000001; //误差

//计算距离
double distance(vector<int>N_1, vector<int>N_2) {
	return double(sqrt(pow((N_1[0] - N_2[0]), 2) + pow((N_1[1] - N_2[1]), 2)));
}
//计算速度
double v_m(int x, double d, double E) {
	return 2 * pow((d-x),3) / (3 * beta*E);
}
double fomugamma_0(int x, double v, double d, double E) {
	return v*E / (d - x) + pow((d - x), 2) / (3 * beta) + H*H / beta;
}
//计算在某一速度下的信息量
double fomu_B(int x, double v, double d, double E) {
	double gamma_0 = fomugamma_0(x, v, d, E);
	double max_B = W / (2 * v)*((d - x)*log2(beta*gamma_0 / (pow((d - x), 2) + H*H)) + alpha*(d - x) / log(2) - alpha*H / log(2)*atan((d - x) / H));
	return max_B;
}
//寻找最优的v
double find_v(double v, double v_max, double B, int x, double d,double E) {
	double v_l = v;
	double v_r = v_max;
	while (v_l < v_r) {
		double mid_v = v_l + (v_r - v_l) / 2;
		double B_op = fomu_B(x, mid_v, d, E);
		if (abs(B_op - B) < err)
			return mid_v;
		else if (B_op < B)
			v_r = mid_v;
		else
			v_l = mid_v;
	}
}
int main()
{
	double B[node_num - 1] = { 1.2*pow(10,6),0.5*pow(10,6),2.1*pow(10,6),1.1*pow(10,6),1.5*pow(10,6),2.1*pow(10,6),1.7*pow(10,6),1.9*pow(10,6),1.3*pow(10,6) };
	double En[node_num - 1] = { 1,1.2,0.7,1.4,1,1.7,2.1,0.3,2.4 };
	srand((unsigned)time(NULL));
	unordered_map<int, vector<int>>nodes;//保存生成各节点的坐标
	vector<vector<double>>dis(node_num,vector<double>(node_num, DBL_MAX));
	vector<vector<double>>t_min(node_num, vector<double>(node_num, DBL_MAX));//保存最短时间矩阵
	nodes[0] = { 0,0 };
	for (int i = 1; i < node_num; ++i) { //随机生成各节点坐标
		nodes[i] = random_coor();
		//cout << nodes[i][1] << endl;
	}
	for (int i = 0; i < node_num; ++i) {
		//cout << "[";
		for (int j = 0; j < node_num; ++j) {
			if (i != j)
				dis[i][j] = distance(nodes[i], nodes[j]);//计算出各节点的距离
			//cout << dis[i][j] << ",";
		}
		//cout << "]" << endl;
	}
	for (int i = 0; i < node_num; ++i) {
		for (int j = 1; j < node_num; ++j) {
			double t_min_m = DBL_MAX;
			for (int x = 0; x <= int(dis[i][j]); ++x) {
				if (i != j) {
					double d = dis[i][j];
					double E = En[j - 1];
					double v_min = v_m(x, d, E);
					double B_max = fomu_B(x, v_min, d, E);
					//cout << B_max << endl;
					if (B_max > B[j - 1]) {
						double v_opt = find_v(v_min, v_max, B[j - 1], x, d, E);
						//cout << v_opt << endl;
						if (int(v_opt)+1>=0) {
							//cout << int(v_opt) << endl;
							double t_op = (d - x) / v_opt + x / v_max;
							//cout << t_op << endl;
							t_min_m = min(t_min_m, t_op);
							//cout << t_min_m << endl;
						}
					}
				}
			}
			t_min[i][j] = t_min_m;
		}
	}
	for (int i = 1; i < node_num; ++i)
		t_min[i][0] = dis[i][0] / v_max;
	for (int i = 0; i < node_num; ++i) {
		cout << "[";
		for (int j = 0; j < node_num; ++j) {
			cout << t_min[i][j] << ",";
		}
		cout << "]"<<endl;
	}
	cout << endl;

	int b = 1 << (node_num - 1);
	double f[1 << (node_num - 1)][node_num];
	int pos[1 << (node_num - 1)][node_num];
	memset(f, DBL_MAX, sizeof(f));
	memset(pos, -1, sizeof(pos));
	for (int i = 0; i < node_num; ++i)
		f[0][i] = t_min[i][0];
	for (int st = 1; st < b - 1; ++st) {
		for (int i = 1; i < node_num; ++i) {
			if (st&(1 << (i - 1))) continue;
			double minn = DBL_MAX;
			for (int j = 1; j < node_num; ++j) {
				if (st&(1 << (j - 1))) {
					double tmp = t_min[i][j] + f[st ^ (1 << (j - 1))][j];
					//cout << tmp << ",";
					if (tmp < minn) {
						minn = tmp;
						f[st][i] = tmp;
						pos[st][i] = j;
					}
				}
			}
		}
	}
	double minn = DBL_MAX;
	for (int k = 1; k < node_num; ++k) {
		double tmp = f[(b - 1) ^ (1 << (k - 1))][k] + t_min[0][k];
		//cout << tmp;
		if (tmp < minn) {
			minn = tmp;
			f[b - 1][0] = tmp;
			pos[b - 1][0] = k;
		}
	}
	//cout << pos[b - 1][0] << endl;
	for (int i = 0; i < node_num; ++i) {
		cout << "[";
		for (int j = 0; j < b; ++j) {
			cout << f[j][i] << ",";
		}
		cout << "]" << endl;
	}
	cout << "res: " << f[b - 1][0];
	system("pause");
    return 0;
}

