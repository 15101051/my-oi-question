#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
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
#define MAXN 5010
struct edge
{
	int to,nxt;
}e[MAXN];
int edgenum = 0;
int lin[MAXN] = {0};
void add(int a,int b)
{
	e[++edgenum] = (edge){b,lin[a]};lin[a] = edgenum;
	return;
}
int ch[MAXN];
int str[MAXN];
int s[MAXN][MAXN];
int dep[MAXN];
int fa[MAXN];
void dfs(int k)
{
	dep[k] = dep[fa[k]] + 1;
	str[dep[k]] = ch[k];
	for(int i = 1;i <= dep[k];++i)s[k][i] = str[i];
	for(int i = lin[k];i != 0;i = e[i].nxt)
	{
		dfs(e[i].to);
	}
	return;
}
#define MOD 998244353
int LCP(int a,int b)
{
	int l = min(dep[a],dep[b]);
	for(int i = 1;i <= l;++i)
	{
		if(s[a][i] != s[b][i])return i - 1;
	}
	return l;
}
int LCS(int a,int b)
{
	int l = min(dep[a],dep[b]);
	for(int i = 1;i <= l;++i)
	{
		if(s[a][dep[a] - i + 1] != s[b][dep[b] - i + 1])return i - 1;
	}
	return l;
}
int main()
{
	freopen("template.in","r",stdin);
	freopen("template.out","w",stdout);
	scanf("%d",&n);
	for(int i = 2;i <= n;++i)add(fa[i] = rd(),i);
	for(int i = 1;i <= n;++i)ch[i] = rd() + 1;
	dfs(1);
	int ans = 0;
	for(int i = 1;i <= n;++i)
	{
		for(int j = i + 1;j <= n;++j)
		{
			ans = (ans + 1ll * LCP(i,j) * LCS(i,j) % MOD) % MOD;
		}
	}
	cout << ans << endl;
	return 0;
}
