#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<map>
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
#define MAXN 100010
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
int fa[MAXN];
int dep[MAXN];
int s[MAXN];
#define MAX 45
#define MAXX 50
#define BASE1 19260817
#define MODN1 1000000007
#define BASE2 1001001
#define MODN2 998244353
typedef unsigned int uint;
int siz[MAXN],son[MAXN];
void dfs1(int k)
{
	siz[k] = 1;
	for(int i = lin[k];i != 0;i = e[i].nxt)
	{
		dfs1(e[i].to);
		siz[k] += siz[e[i].to];
		if(son[k] == 0 || siz[e[i].to] > siz[son[k]])son[k] = e[i].to;
	}
	return;
}
#define pii pair<uint,uint>
#define fi first
#define se second
#define mp make_pair
long long kkk = 0;
pii hash[MAXN][MAXX];
typedef unsigned long long ull;
ull hsh[MAXN][MAXX];
struct table
{
	#define MO 1001001
	int head[MO];
	ull st[MAXN];
	int val[MAXN],nxt[MAXN];
	int cntnum;
	int stack[MAXN],top;
	int& operator [] (ull x)
	{
		int modx = x % MO;
		for(int i = head[modx];i != 0;i = nxt[i])
		{
			if(st[i] == x)return val[i];
		}
		stack[++top] = modx;
		++cntnum;
		st[cntnum] = x;
		val[cntnum] = 0;
		nxt[cntnum] = head[modx];
		head[modx] = cntnum;
		return val[cntnum];
	}
	void clear()
	{
		while(top)head[stack[top--]] = 0;
		cntnum = 0;
		return;
	}
};
namespace T
{
	table f[MAXX];
	table sum[MAXX];
	int ans = 0;
	void clear()
	{
		for(int i = 0;i < MAX;++i)f[i].clear();
		for(int i = 0;i < MAX;++i)sum[i].clear();
		ans = 0;
		return;
	}
	void query(int k)
	{//cout << "query " << k << endl;
		for(int i = 0;i <= min(MAX - 1,dep[k]);++i)
		{//cout << i << " : " << hash[k][i] << " " << hash[k][i + 1] << " " << (sum[i][hash[k][i]] - f[i][mp(hash[k][i],hash[k][i + 1])]) * i << endl;
			ans += (sum[i][hsh[k][i]] - f[i][hsh[k][i] ^ hsh[k][i + 1]]) * i;
		}
		return;
	}
	void insert(int k)
	{//cout << "insert " << k << endl;
		for(int i = 0;i <= min(MAX - 1,dep[k]);++i)
		{//cout << i << " : " << hash[k][i] << " " << hash[k][i + 1] << endl;
			++sum[i][hsh[k][i]];++kkk;
			++f[i][hsh[k][i] ^ hsh[k][i + 1]];
		}
		return;
	}
}
void calc(int k)
{
	T::query(k);
	for(int i = lin[k];i != 0;i = e[i].nxt)calc(e[i].to);
	return;
}
void add(int k)
{
	T::insert(k);
	for(int i = lin[k];i != 0;i = e[i].nxt)add(e[i].to);
	return;
}
int ans = 0;
#define MOD 998244353
void dfs(int k,int ty)
{
	if(son[k] != 0)
	{
		for(int i = lin[k];i != 0;i = e[i].nxt)if(e[i].to != son[k])dfs(e[i].to,0);
		dfs(son[k],1);
		T::ans = 0;
		for(int i = lin[k];i != 0;i = e[i].nxt)
		{
			if(e[i].to == son[k])continue;
			calc(e[i].to);
			add(e[i].to);
		}
	}
	T::query(k);
	T::insert(k);
	ans = (ans + 1ll * T::ans * dep[k] % MOD) % MOD;
	if(ty == 0)T::clear();
	return;
}
int main()
{
	freopen("template.in","r",stdin);
	freopen("template.out","w",stdout);
	scanf("%d",&n);
	for(int i = 2;i <= n;++i)add(fa[i] = rd(),i);
	for(int i = 1;i <= n;++i)dep[i] = dep[fa[i]] + 1;
	for(int i = 1;i <= n;++i)s[i] = rd() + 1;
	s[0] = 257;
	int cnt = -1;
	for(int i = 1;i <= n;++i)
	{
		for(int j = 1,cur = i;j <= min(dep[i],MAX);++j,cur = fa[cur])
		{
			hash[i][j].fi = (1ll * hash[i][j - 1].fi * BASE1 + s[cur]) % MODN1;
			hash[i][j].se = (1ll * hash[i][j - 1].se * BASE2 + s[cur]) % MODN2;
			hsh[i][j] = (1ull * hash[i][j].fi) << 32 | hash[i][j].se;
		}
		for(int j = dep[i] + 1;j <= MAX;++j)hsh[i][j] = --cnt;
	}
	dfs1(1);
	dfs(1,0);
	cout << ans << endl;
	fclose(stdin);
	fclose(stdout);
	return 0;
}
