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
#define MAXN 500010
#define LOG 19
#define MAXM 256
int fa[MAXN],ch[MAXN];
#define R register
#define I inline
struct edge
{
	int to,nxt;
}e[MAXN];
int edgenum = 0;
int lin[MAXN] = {0};
I void add(int a,int b)
{
	e[++edgenum] = (edge){b,lin[a]};lin[a] = edgenum;
	return;
}
int dep[MAXN]; 
int f[MAXN][LOG];
int sa[MAXN],rnk[LOG][MAXN],h[MAXN],c[MAXN],c1[MAXN],c2[MAXN],rk[MAXN];
vector<int> v[MAXN];
struct table
{
	int head[MAXN],nxt[MAXN],val[MAXN];
	int cntnum;
	void add(int a,int b)
	{
		++cntnum;val[cntnum] = b;nxt[cntnum] = head[a];head[a] = cntnum;
		return;
	}
	table(){cntnum = 0;}
}l;
I void make_SA(int n,int m)
{
	R int p;
	R int *x = c1,*y = c2;
	for(R int i = 1;i <= m;++i)c[i] = 0;
	for(R int i = 1;i <= n;++i)++c[x[i] = ch[i]];
	for(R int i = 1;i <= m;++i)c[i] += c[i - 1];
	for(R int i = n;i >= 1;--i)sa[c[x[i]]--] = i;
	for(R int i = 1;i <= n;++i)rnk[0][i] = ch[i];
	for(R int k = 1,j = 1;(1 << j) <= n;k = k << 1,++j)
	{
		p = 0;
		l.cntnum = 0;
		for(R int i = 1;i <= n;++i)l.head[i] = 0;
		for(R int i = 1;i <= m;++i)c[i] = 0;
		for(R int i = 1;i <= n;++i)if(dep[sa[i]] > k)l.add(f[sa[i]][j - 1],sa[i]);
		for(R int i = 1;i <= n;++i)for(R int k = l.head[sa[i]];k != 0;k = l.nxt[k])y[++p] = l.val[k];
		for(R int i = 1;i <= n;++i)if(dep[i] <= k)y[++p] = i;
		for(R int i = 1;i <= n;++i)++c[x[y[i]]];
		for(R int i = 1;i <= m;++i)c[i] += c[i - 1];
		for(R int i = n;i >= 1;--i)sa[c[x[y[i]]]--] = y[i];
		p = 1;
		swap(x,y);
		x[sa[1]] = 1;
		for(R int i = 2;i <= n;++i)
		{
			x[sa[i]] = (y[sa[i]] == y[sa[i - 1]] && y[f[sa[i]][j - 1]] == y[f[sa[i - 1]][j - 1]] ? p : ++p);
		}
		for(R int i = 1;i <= n;++i)rnk[j][i] = x[i];
		m = p;
	}
	for(R int i = 1;i <= n;++i)rk[sa[i]] = i;
	return;
}
I int LCP(int a,int b)
{
	if(a == b)return dep[a];
	R int res = 0;
	for(R int i = LOG - 1;i >= 0;--i)
	{
		if(min(dep[a],dep[b]) >= (1 << i) && rnk[i][a] == rnk[i][b])
		{
			res += (1 << i);
			a = f[a][i];b = f[b][i];
		}
	}
	return res;
}
#define pii pair<int,int>
#define fi first
#define se second
pii st[MAXN][LOG];
int lg[MAXN];
I pii mymin(pii a,pii b)
{
	if(a.fi != b.fi)return min(a,b);
	else if(rand() % 2 == 0)return a;
	else return b;
}
I pii query(int l,int r)
{
	R int len = lg[r - l + 1];
	return mymin(st[l][len],st[r - (1 << len) + 1][len]);
}
int trt,lc[MAXN << 1],rc[MAXN << 1],len[MAXN << 1],tot;
long long mv[MAXN];
int le[MAXN];
void build(int &rt,int l,int r)
{
	if(l == r)
	{
		rt = l;
		return;
	}
	rt = ++tot;
	R int pos = query(l + 1,r).second;
	len[rt] = query(l + 1,r).first;
	for(R int i = l;i < pos;++i)mv[i] = mv[i] << 1 | 0,++le[i];
	build(lc[rt],l,pos - 1);
	for(R int i = pos;i <= r;++i)mv[i] = mv[i] << 1 | 1,++le[i];
	build(rc[rt],pos,r);
	return;
}
struct node
{
	int lc,rc;
	int sum,l;
	node(){sum = 0;}
}t[MAXN * 30];
int ptr = 0;
I int newnode(){return ++ptr;}
int root[MAXN];
I void insert(int k)
{
	root[k] = newnode();
	t[root[k]].sum = 1;t[root[k]].l = len[trt];
	R int cur = root[k],cur_ = trt;
	k = rk[k];
	for(R int i = le[k] - 1;i >= 0;--i)
	{
		if((mv[k] >> i) & 1)
		{
			cur = t[cur].rc = newnode();
			cur_ = rc[cur_];
		}
		else
		{
			cur = t[cur].lc = newnode();
			cur_ = lc[cur_];
		}
		t[cur].sum = 1;
		t[cur].l = len[cur_];
	}
	return;
}
#define MOD 998244353
int ans = 0;
int merge(int x,int y,int val)
{
	if(x == 0 || y == 0)return x + y;
	ans = (ans + 1ll * val * t[x].l * t[t[x].lc].sum % MOD * t[t[y].rc].sum % MOD) % MOD;
	ans = (ans + 1ll * val * t[x].l * t[t[x].rc].sum % MOD * t[t[y].lc].sum % MOD) % MOD;
	t[x].sum += t[y].sum;
	t[x].lc = merge(t[x].lc,t[y].lc,val);
	t[x].rc = merge(t[x].rc,t[y].rc,val);
	return x;
}
void dfs(int k)
{
	for(R int i = lin[k];i != 0;i = e[i].nxt)
	{
		dfs(e[i].to);
		root[k] = merge(root[k],root[e[i].to],dep[k]);
	}
	return;
}
int main()
{
	freopen("template.in","r",stdin);
	freopen("template.out","w",stdout);
	scanf("%d",&n);
	for(R int i = 2;i <= n;++i)add(fa[i] = rd(),i),lg[i] = lg[i >> 1] + 1;
	for(R int i = 1;i <= n;++i)dep[i] = dep[fa[i]] + 1,ch[i] = rd() + 1;
	for(R int i = 1;i <= n;++i)f[i][0] = fa[i];
	for(R int k = 1;k < LOG;++k)
		for(R int i = 1;i <= n;++i)f[i][k] = f[f[i][k - 1]][k - 1];
	make_SA(n,MAXM);
	for(R int i = 1;i <= n;++i)st[i][0] = make_pair(h[i] = LCP(sa[i],sa[i - 1]),i);
	for(R int k = 1,l = 1;k < LOG;++k,l = l << 1)
		for(R int i = 1;i <= n - (l << 1) + 1;++i)st[i][k] = mymin(st[i][k - 1],st[i + l][k - 1]);
	tot = n;
	build(trt,1,n);
	for(R int i = 1;i <= n;++i)insert(i);
	dfs(1);
	printf("%d\n",ans);
	fclose(stdin);
	fclose(stdout);
	return 0;
}
