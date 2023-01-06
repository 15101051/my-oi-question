#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
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
#define LOG 18
#define MAXM 256
int fa[MAXN],ch[MAXN];
int dep[MAXN]; 
int f[MAXN][LOG];
int sa[MAXN],rnk[LOG][MAXN],h[MAXN],c[MAXN],c1[MAXN],c2[MAXN];
int rk[MAXN];
vector<int> v[MAXN];
void make_SA(int n,int m)
{
	int p;
	int *x = c1,*y = c2;
	for(int i = 1;i <= m;++i)c[i] = 0;
	for(int i = 1;i <= n;++i)++c[x[i] = ch[i]];
	for(int i = 1;i <= m;++i)c[i] += c[i - 1];
	for(int i = n;i >= 1;--i)sa[c[x[i]]--] = i;
	for(int i = 1;i <= n;++i)rnk[0][i] = ch[i];
	for(int k = 1,j = 1;(1 << j) <= n;k = k << 1,++j)
	{
		p = 0;
		for(int i = 1;i <= n;++i)v[i].clear();
		for(int i = 1;i <= m;++i)c[i] = 0;
		for(int i = 1;i <= n;++i)if(dep[sa[i]] > k)v[f[sa[i]][j - 1]].push_back(sa[i]);
		for(int i = 1;i <= n;++i)for(vector<int>::iterator it = v[sa[i]].begin();it != v[sa[i]].end();++it)y[++p] = *it;
		for(int i = 1;i <= n;++i)if(dep[i] <= k)y[++p] = i;
		for(int i = 1;i <= n;++i)++c[x[y[i]]];
		for(int i = 1;i <= m;++i)c[i] += c[i - 1];
		for(int i = n;i >= 1;--i)sa[c[x[y[i]]]--] = y[i];
		p = 1;
		swap(x,y);
		x[sa[1]] = 1;
		for(int i = 2;i <= n;++i)
		{
			x[sa[i]] = (y[sa[i]] == y[sa[i - 1]] && y[f[sa[i]][j - 1]] == y[f[sa[i - 1]][j - 1]] ? p : ++p);
		}
		for(int i = 1;i <= n;++i)rnk[j][i] = x[i];
		m = p;
	}
	for(int i = 1;i <= n;++i)rk[sa[i]] = i;
	return;
}
int LCP(int a,int b)
{
	if(a == b)return dep[a];
	int res = 0;
	for(int i = LOG - 1;i >= 0;--i)
	{
		if(min(dep[a],dep[b]) >= (1 << i) && rnk[i][a] == rnk[i][b])
		{
			res += (1 << i);
			a = f[a][i];b = f[b][i];
		}
	}
	return res;
}
struct ST
{
	pair<int,int> f[MAXN << 2][LOG];
	int lg[MAXN << 1];
	int s;
	void build(int n)
	{
		for(int i = 2;i <= n;++i)lg[i] = lg[i >> 1] + 1;
		s = n;
		for(int k = 1,l = 1;k < LOG;++k,l = l << 1)
			for(int i = 1;i <= s - (l << 1) + 1;++i)f[i][k] = min(f[i][k - 1],f[i + l][k - 1]);
		return;
	}
	pair<int,int> query(int l,int r)
	{
		if(l > r)swap(l,r);
		int len = lg[r - l + 1];
		return min(f[l][len],f[r - (1 << len) + 1][len]);
	}
}lca,hei;
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
int dfn[MAXN << 1],pos[MAXN],tot = 0;
void dfs(int k)
{
	dfn[++tot] = k;pos[k] = tot;
	for(int i = lin[k];i != 0;i = e[i].nxt)
	{
		dfs(e[i].to);
		dfn[++tot] = k;
	}
	return;
}
#define MOD 998244353
int lcp(int a,int b)
{
	return dep[lca.query(pos[a],pos[b]).second];
}
int lcs(int a,int b)
{
	if(rk[a] > rk[b])swap(a,b);
	return hei.query(rk[a] + 1,rk[b]).first;
}
int main()
{
	freopen("template.in","r",stdin);
	freopen("template.out","w",stdout);
	scanf("%d",&n);
	for(int i = 2;i <= n;++i)add(fa[i] = rd(),i);
	for(int i = 1;i <= n;++i)dep[i] = dep[fa[i]] + 1;
	for(int i = 1;i <= n;++i)ch[i] = rd() + 1;
	for(int i = 1;i <= n;++i)f[i][0] = fa[i];
	for(int k = 1;k < LOG;++k)
		for(int i = 1;i <= n;++i)f[i][k] = f[f[i][k - 1]][k - 1];
	make_SA(n,MAXM);
	for(int i = 1;i <= n;++i)hei.f[i][0].first = h[i] = LCP(sa[i],sa[i - 1]);
	hei.build(n);
	dfs(1);
	for(int i = 1;i <= tot;++i)lca.f[i][0] = make_pair(dep[dfn[i]],dfn[i]);
	lca.build(tot);
	int ans = 0;
	for(int i = 1;i <= n;++i)
	{
		for(int j = i + 1;j <= n;++j)
		{
			ans = (ans + 1ll * lcs(i,j) * lcp(i,j) % MOD) % MOD;
		}
	}
	cout << ans << endl;
	fclose(stdin);
	fclose(stdout);
	return 0;
}
