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
#define MAXN 300010
int nxt[MAXN];
int s[MAXN];
int dfn[MAXN];
int sa[MAXN],rnk[MAXN],c[MAXN],c1[MAXN],c2[MAXN],h[MAXN];
void make_SA(int n,int m)
{
	int p = 0;
	int *x = c1,*y = c2;
	for(int i = 1;i <= m;++i)c[i] = 0;
	for(int i = 1;i <= n;++i)++c[x[i] = s[i]];
	for(int i = 1;i <= m;++i)c[i] += c[i - 1];
	for(int i = n;i >= 1;--i)sa[c[x[i]]--] = i;
	for(int k = 1;k <= n;k = k << 1)
	{
		p = 0;
		for(int i = 1;i <= n;++i)y[i] = 0;
		for(int i = 1;i <= m;++i)c[i] = 0;
		for(int i = n - k + 1;i <= n;++i)y[++p] = i;
		for(int i = 1;i <= n;++i)if(sa[i] > k)y[++p] = sa[i] - k;
		for(int i = 1;i <= n;++i)++c[x[y[i]]];
		for(int i = 1;i <= m;++i)c[i] += c[i - 1];
		for(int i = n;i >= 1;--i)sa[c[x[y[i]]]--] = y[i];
		p = 1;
		swap(x,y);
		x[sa[1]] = 1;
		for(int i = 2;i <= n;++i)
		{
			x[sa[i]] = (y[sa[i]] == y[sa[i - 1]] && y[sa[i] + k] == y[sa[i - 1] + k] ? p : ++p);
		}
		if(p >= n)break;
		m = p;
	}
	for(int i = 1;i <= n;++i)rnk[sa[i]] = i;
	int k = 0;
	for(int i = 1;i <= n;++i)
	{
		if(rnk[i] == 1)continue;
		if(k)--k;
		int j = sa[rnk[i] - 1];
		while(j + k <= n && i + k <= n && s[j + k] == s[i + k])++k;
		h[rnk[i]] = k;
	}
	return;
}
int root,lc[MAXN << 1],rc[MAXN << 1],fa[MAXN << 1],len[MAXN << 1];
int sum[MAXN << 1];
int tot;
int build(int &rt,int l,int r)
{
	if(l == r)
	{
		rt = l;
		return rt;
	}
	int mn = 0x3f3f3f3f,pos = 0;
	for(int i = l + 1;i <= r;++i)if(h[i] < mn){pos = i;mn = h[i];}
	rt = ++tot;len[rt] = h[pos];
	fa[build(lc[rt],l,pos - 1)] = rt;
	fa[build(rc[rt],pos,r)] = rt;
	return rt;
}
#define MOD 998244353
int main()
{
	freopen("template.in","r",stdin);
	freopen("template.out","w",stdout);
	scanf("%d",&n);
	for(int i = 2;i <= n;++i)nxt[rd()] = i;
	for(int i = 1,cur = 1;i <= n;++i)
	{
		dfn[cur] = i;
		cur = nxt[cur];
	}
	for(int i = 1;i <= n;++i)s[n - dfn[i] + 1] = rd() + 1;
	make_SA(n,256);
	tot = n;
	build(root,1,n);
	int ans = 0;
	for(int i = 2;i <= n;++i)
	{
		int cur = rnk[i - 1];
		while(cur != 0)
		{
			++sum[cur];
			cur = fa[cur];
		}
		cur = rnk[i];
		while(fa[cur] != 0)
		{
			ans = (ans + 1ll * len[fa[cur]] * (n - i + 1) % MOD * (sum[fa[cur]] - sum[cur]) % MOD) % MOD;
			cur = fa[cur];
		}
	}
	cout << ans << endl;
	fclose(stdin);
	fclose(stdout);
	return 0;
}
