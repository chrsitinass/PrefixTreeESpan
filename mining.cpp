#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

#define make_pair MP

typedef vector<int> Tree;

int total = 0;
// double thres[] = {3.5, 3, 2.5, 2, 1.75, 1.7, 1.65, 1.6, 1.55, 1.5};
double thres[] = {1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
int total_pattern = 0;
int cur_min;
int MIN_SUPPORT[15];

struct Element{
	int label, attach;
	Element(int l, int a): label(l), attach(a) {}
	bool operator < (const Element & x) const {
		return label < x.label;
	}
	void print() const {
		printf("Element(label, attach): %d %d\n", label, attach);
	}
};
struct Proj {
	vector<int> position;
	int tree_id;
	Proj(int id, int pos) {
		tree_id = id;
		position = {pos};
	}
};
vector<Tree> data;	// database
vector<vector<int>> end_flag;	// end flag mapping
unordered_map<int, int> label_count;

/* 
 *	Get test data
*/
void init() {
	string path[] = {
		"./treedata/small.data",
		"./treedata/CSlog.data", 
		"./treedata/D10.data", 
		"./treedata/F5.data", 
		"./treedata/T1M.data"
	};
	freopen(path[1].c_str(), "r", stdin);
	int t = 0, sum = 0;
	vector<int> cur;
	set<int> occur_in_tree;
	int total_num = 0;
	while (scanf("%d", &t) != EOF) {
		cur.push_back(t);
		total_num++;
		if (t < 0) sum -= 1;
		else {
			sum += 1;
			occur_in_tree.insert(t);
		}
		if (sum == 0) {
			data.push_back(cur);
			cur.clear();
			++total;
			for (auto it: occur_in_tree) {
				label_count[it] += 1;
			}
			occur_in_tree.clear();
		}
	}
	printf("The size of database: %d\n", total);
	printf("The average size of trees in database: %lf\n", total_num * 1.0 / total);
	for (int i = 0; i < 10; i++) {
		MIN_SUPPORT[i] = thres[i] * total / 100;
	}
}

/*
 * Find each node's end flag '-1'
*/
void getTreeStruct() {
	for (auto tree: data) {
		stack<int> sta;
		vector<int> temp;
		for (auto i = 0; i < tree.size(); i++) {
			if (tree[i] >= 0) {
				sta.push(i);
				temp.push_back(-1);
			} else {
				temp.push_back(sta.top());
				sta.pop();
			}
		}
		end_flag.push_back(temp);
	}
}

/*
 * Scan ProDB(D,S) once to find all frequent GEs.
*/
vector<Element> getGrowthElement(const vector<Proj> &projection) {
	vector<Element> growth_element;
	// occurrences of candidate growth element in projected database
	map<Element, set<int>> candidate_count;
	for (auto it: projection) {
		int id = it.tree_id;
		vector<int> pos = it.position;
		int current_attach = pos.size() - 1;
		for (auto i = pos[current_attach] + 1; i < data[id].size(); i++) {
			int t = data[id][i];
			if (t == -1) {
				while (current_attach >= 0 && end_flag[id][i] <= pos[current_attach]) {
					--current_attach;
				}
				if (current_attach < 0) break;
			} else {
				candidate_count[Element(t, current_attach)].insert(id);
			}
		}
	}
	for (auto it: candidate_count) {
		if (it.second.size() >= MIN_SUPPORT[cur_min]) {
			growth_element.push_back(it.first);
			// it.first.print();
		}
	}
	return growth_element;
}

/*
 * Extent pattern by ge to form a subtree pattern `new_pattern`, 
 * and output `new_pattern` 
*/
Tree generatePattern(Tree pattern, Element ge) {
	stack<int> sta;
	Tree res;
	for (auto i = 0; i < pattern.size(); i++) {
		int t = pattern[i];
		if (t == -1) {
			if (sta.top() == ge.attach) {
				res.push_back(ge.label);	// attach growth element
				res.push_back(-1);
			}
			res.push_back(-1);	// end flag of the origin subtree
			sta.pop();
		} else {
			sta.push(i);
			res.push_back(t);
		}
	}
	return res;
}

/*
 * Find all Occurrences of ge in ProDB(D,pattern), 
 * and construct <new_pattern> -projected database
 * through collecting all corresponding Project Instances in ProDB(D,S)
*/
vector<Proj> generateProj(Element ge, const vector<Proj> &projection) {
	vector<Proj> res;
	for (auto it: projection) {
		int id = it.tree_id;
		vector<int> &pos = it.position;
		int current_attach = pos.size() - 1;
		for (auto i = pos[current_attach] + 1; i < data[id].size(); i++) {
			int t = data[id][i];
			if (t == ge.label && current_attach == ge.attach) {
				auto new_it = it;
				new_it.position.push_back(i);
				res.push_back(new_it);
			}
			if (t == -1) {
				while (current_attach >= 0 && end_flag[id][i] <= pos[current_attach]) {
					--current_attach;
				}
				if (current_attach < 0) break;
			}
		}
	}
	return res;
}

void fre(Tree pattern, const vector<Proj> &projection) {
	if (projection.size() == 0) {
		return;
	}
	vector<Element> growth_element = getGrowthElement(projection);
	for (auto ge: growth_element) {
		Tree new_pattern = generatePattern(pattern, ge);
		// Output new_pattern into a file as a result
		const vector<Proj> &new_projection = generateProj(ge, projection);
		fre(new_pattern, new_projection);
	}
}

/*
 * We don't need to record every subtree pattern
 */
void fre(const vector<Proj> &projection) {
	if (projection.size() == 0) {
		return;
	}
	const vector<Element> &growth_element = getGrowthElement(projection);
	total_pattern += growth_element.size();
	for (auto ge: growth_element) {
		const vector<Proj> &new_projection = generateProj(ge, projection);
		fre(new_projection);
	}
}

void prefixTreeESpan() {
	for (auto it: label_count) {
		if (it.second >= MIN_SUPPORT[cur_min]) {
			total_pattern += 1;
			// start with frequent label
			int label = it.first;
			Tree pattern = {label, -1};
			// init <label -1> projection
			vector<Proj> projection;
			for (auto id = 0; id < data.size(); id++) {
				for (auto pos = 0; pos < data[id].size(); pos++) {
					int num = data[id][pos];
					if (num == label) {
						projection.push_back(Proj(id, pos));
					}
				}
			}
			// mining <label -1>
			fre(pattern, projection);
			// fre(projection);
		}
	}
}

int main() {
	init();
	getTreeStruct();
	for (cur_min = 0; cur_min < 10; cur_min++) {
		total_pattern = 0;
		clock_t start = clock();
		prefixTreeESpan();
		clock_t end = clock();
		printf("Mining finished.\nTotal time cost: %lf seconds.\nTotal pattern found: %d\n", 
				(double)(end - start) / CLOCKS_PER_SEC, total_pattern);
	}
	return 0;
}