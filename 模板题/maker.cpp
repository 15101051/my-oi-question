#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<map>
#include<ctime>
#include<cctype>
#include<cstring>
using namespace std;
inline int rd()
{
	register int res = 0,f = 1;register char c = getchar();
	while(!isdigit(c)){if(c == '-')f = -1;c = getchar();}
	while(isdigit(c))res = (res << 1) + (res << 3) + c - '0',c = getchar();
	return res * f;
}
int n;
int sigma = 2;
map<int,int> check[500010];
int anc[500010];
int randchar(int i)
{
    int ch = rand() * rand() % sigma;
    while(check[anc[i]][ch])ch = rand() * rand() % sigma;
    check[anc[i]][ch] = 1;
    return ch;
}
int deg[500010];
int getfa(int i)
{
	if(i <= 20)
	{
		int fa = rand() % (i - 1) + 1;
		while(deg[fa] == 2)fa = rand() % (i - 1) + 1;
		++deg[fa];
		return fa;
	}
	else
	{
		int fa = i - rand() % 20 - 1;
		while(deg[fa] == 2)fa = i - rand() % 20 - 1;
		++deg[fa];
		return fa;
	}
}
int main()
{
	freopen("template10.in","w",stdout);
	srand(time(NULL));
	n = 50000;
	cout << n << endl;
	//for(int i = 2;i <= n;++i)cout << (anc[i] = i - 1) << " ";
	for(int i = 2;i <= n;++i)
	{
		cout << (anc[i] = getfa(i)) << " ";
	}
	cout << endl;
	for(int i = 1;i <= n;++i)cout << randchar(i) << " ";cout << endl;
	return 0;
}
