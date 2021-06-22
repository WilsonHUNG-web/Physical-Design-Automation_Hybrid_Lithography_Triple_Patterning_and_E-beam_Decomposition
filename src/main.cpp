#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>  
#include <chrono>
#include <map>
#include <list>
#include <set>

#define MAX_COORD 10000000;
#define delta 0.000001; //for dummy extension
#define infinity 100000000; //for infinite distance in Dijk

using namespace std;

//global var 
double dmin;
int ccccc = 0;
int num_poly = 0;
int numAs = 0; //number of isolated polygons

int c1=0;
int c2=0;
int c3=0;
int c4=0;

int TotalnumNode = 0;
int TotalnumEdge = 0;
//

//class coord {//no use
//public:
//	double x;
//	double y;
//	//constructor
//	coord(double X, double Y) {
//		this->x = X;
//		this->y = Y;
//	}
//
//};

//for data in/out
class rect {
public:
	double xl, yl, xr, yr;
	double e_xl, e_yl, e_xr, e_yr; //extended size with dmin


	bool is_rect_conflict(rect* RE) {//check extended this->rect is overlapping with another rect RE
		if (((RE->xl < this->e_xl && RE->xr < this->e_xl) || (RE->xl > this->e_xr && RE->xr > this->e_xr)) ||
			((RE->yl < this->e_yl && RE->yr < this->e_yl) || (RE->yl > this->e_yr && RE->yr > this->e_yr)))
		{
			//cout << "false" << endl;
			return false;
		}
		//if ((RE->xl < this->e_xl && RE->xr < this->e_xl)) cout << "ffadfsdfd" << endl;
		return true;
	}
	//constructor
	rect(double xl_in, double yl_in, double xr_in, double yr_in) {
		this->xl = xl_in;
		this->xr = xr_in;
		this->yl = yl_in;
		this->yr = yr_in;

		this->e_xl = xl_in - dmin;
		this->e_yl = yl_in - dmin;
		this->e_xr = xr_in + dmin;
		this->e_yr = yr_in + dmin;
	}
};

class poly {
public:
	string name;
	double max_x, max_y; //right-up boundary
	double min_x, min_y; //left-down boundary
	double dummy_x; //dummy extension coordinate
	//vector<coord> outline;
	char color;
	poly* ptr;
	vector<poly*> CG; //constraint graph
	vector<rect*> rects;
	bool inSG;


	void set_to_A() {
		this->color = 'A';
	}
	void set_to_B() {
		this->color = 'B';
	}
	void set_to_C() {
		this->color = 'C';
	}

	void set_dummy_ext() {
		double rightmost_x = 0;
		bool flag = false;
		for (unsigned int i = 0; i < this->CG.size(); i++) {
			double cur_min_x = this->CG[i]->min_x;
			if (cur_min_x > this->min_x) { //if another poly on CG is on the right to this poly
				if (cur_min_x > rightmost_x) {
					rightmost_x = cur_min_x;
					flag = true;
				}
			}
		}
		if (flag) { this->dummy_x = rightmost_x - delta; }
		else { this->dummy_x = this->max_x; }
	}

	void push_rect(rect* RE) {//used
		if (RE->xr > this->max_x) this->max_x = RE->xr;
		if (RE->yr > this->max_y) this->max_y = RE->yr;
		if (RE->xl < this->min_x) this->min_x = RE->xl;
		if (RE->yl < this->min_y) this->min_y = RE->yl;

		this->dummy_x = this->max_x;
		rects.push_back(RE);
		//this->outline.push_back(coordinate);
	}

	bool is_in_rect_boundary(poly* p) {
		bool x_ok = false, y_ok = false;
		//check x direction
		if (this->max_x + dmin< p->min_x - dmin || this->min_x - dmin > p->max_x + dmin) x_ok = true;
		//check y direction 
		if (this->max_y + dmin< p->min_y - dmin || this->min_y - dmin > p->max_y + dmin) y_ok = true;

		if (x_ok || y_ok) return true;
		else return false;
	}

	bool is_in_cg(poly* p) {
		for (size_t i = 0; i < this->rects.size(); i++) { //i for this
			for (size_t j = 0; j < p->rects.size(); j++) { //j for p
				if (this->rects[i]->is_rect_conflict(p->rects[j])) {
					this->addToCG(p);
					//cout << "CG size: " << this->CG.size() << endl;
					return true;
				}
			}

		}

		return false;
	}


	void addToCG(poly* conf_poly) {
		this->CG.push_back(conf_poly);
	}


	//constructor
	poly(char color, string name) {
		this->name = name;
		this->max_x = this->max_y = 0;
		this->min_x = this->min_y = MAX_COORD;
		this->color = color;
		this->inSG = true;
	}



};

//for graph

class Element {
public:
	poly* polygon;
	int color; //1 to 4 for A to E
	Element(poly* p, int c) : polygon(p), color(c) {}

};

class Node {
public:
	vector<Element*> elements;
	int idx_adjList;

	int get_num_E() {
		int count_E = 0;
		for (unsigned int i = 0; i < elements.size(); i++) {
			if (this->elements[i]->color == 'E') count_E++;
		}
		return count_E;
	}

	Node(int idx): idx_adjList(idx){}
};

class adjNode {
public:
	Node* key;
	Node* next;
	int cost;
	adjNode(Node* src, Node* des, int E_num_in_des): key(src), next(des), cost(E_num_in_des){}
};

class NodeGroup {
public:
	vector<poly*> polygon_set;
	vector<Node*> nodes;

	NodeGroup(vector<poly*> set_from_cutline_set) {
		this->polygon_set = set_from_cutline_set;
	}
};

class Edge {
public:
	Node* src;
	Node* des;
	int cost;
	Edge(int num_of_E_in_des) : cost(num_of_E_in_des) {}

};
//


class row {
public:
	string name;
	double x, y;
	vector<poly*> polygons;
	vector<poly*> sorted_poly;
	vector<vector<poly*>> cutline_def;//put the poly with ths same min_x to the same vector<poly*>
	//each vec in this 'cutline_def' is an s-star set 
	vector<vector<poly*>> cutline_set;
	int node_num;

	//for graph
	vector<NodeGroup*> nodegroups;
	double get_polygon_dist(poly p1, poly p2);

	//setup the sorted vector and delete no-cg poly from it
	void sort_polygons() { //ascending order
		sort(this->sorted_poly.begin(), this->sorted_poly.end(), [](poly* const &A, poly* const&B) { return A->min_x < B->min_x; });
	}

	//bool poly_is_not_in_cg(poly* p) {//for remove_if()
	//	return p->inSG;
	//}

	//void del_noCG_poly() {
	//	vector<poly*>::iterator it;
	//	it = remove_if(this->sorted_poly.begin(), this->sorted_poly.end(), this->poly_is_not_in_cg);
	//}
	//

	vector<poly*> copy_vec_poly_in_cg(vector<poly*> ps) {
		vector<poly*> vec_new;
		for (size_t i = 0; i < ps.size(); i++) {
			if (ps[i]->inSG) vec_new.push_back(ps[i]);
		}
		return vec_new;
	}

	void set_dummy_ext_row() {
		for (unsigned int i = 0; i < this->sorted_poly.size(); i++) {
			this->sorted_poly[i]->set_dummy_ext();
		}
	}

	void set_cutline() {
		//cout << "g0 ";
		this->sorted_poly = this->copy_vec_poly_in_cg(this->polygons);
		//cout << "g1 ";
		this->sort_polygons();
		//cout << "g2 ";
		double cur_x = -1, prev_x = -1;
		//cout << this->name << endl;
		for (unsigned int i = 0; i < this->sorted_poly.size(); i++) {//put the poly with ths same min_x to the same vector
			cur_x = sorted_poly[i]->min_x;
			if (cur_x == prev_x) {
				this->cutline_def.back().push_back(this->sorted_poly[i]);
				//cout << ++ccccc << endl;
			}
			else {
				vector<poly*> vec;
				vec.push_back(this->sorted_poly[i]);
				this->cutline_def.push_back(vec);
			}

			prev_x = sorted_poly[i]->min_x;
			//cutline_def.set_dummy_ext()
		}

	}

	void create_cutline_set() {
		this->cutline_set = this->cutline_def;

		for (unsigned int i = 1; i < this->cutline_set.size(); i++) { //dont need to do i =0
			double cut_x = this->cutline_set[i][0]->min_x;
			for (unsigned int j = 0; j < i; j++) {
				for (unsigned int k = 0; k < cutline_def[j].size(); k++) {  //use original cutline instead of cutline_set
					poly* temp = cutline_set[j][k];
					if (temp->dummy_x > cut_x) cutline_set[i].push_back(temp);
				}
			}
		}

	}

	void build_cg() {
		for (unsigned int i = 0; i < this->polygons.size(); i++) {

			for (unsigned int j = 0; j < this->polygons.size(); j++) {
				if (j == i) continue;
				this->polygons[i]->is_in_cg(this->polygons[j]);
			}
			if (this->polygons[i]->CG.size() == 0) { //if poly not in any CG, assign it with color A and don't add to SG
				this->polygons[i]->set_to_A();
				this->polygons[i]->inSG = false;
				numAs++;
			}
		}
		//cout  <<endl;
	}

	void create_elements_n_nodes() {
		for(unsigned int i = 0;i < this->cutline_set.size(); i++)
			if (this->cutline_set[i].size()==1 ){
				for (int j = 0; j < 4; j++) {


			}
			}
			else if (this->cutline_set[i].size() == 2){
				for (int j = 0; j < 2; j++) {

				}
			}
			else if (this->cutline_set[i].size() == 3){}
			else if (this->cutline_set[i].size() == 4){}

			
				
			

	}

	void create_nodes() {
		for (unsigned int i = 0; i < this->cutline_set.size(); i++) {
			NodeGroup* temp_ng = new NodeGroup(cutline_set[i]);
			int numPoly_in_set = cutline_set[i].size();
			if(numPoly_in_set == 1) {
				c1++;
				for (int j = 1; j < 5; j++) {
					Node* temp_n = new Node(this->node_num++);
					//cout << temp_n->idx_adjList << endl;
					poly* temp_p = cutline_set[i][0];
					Element* temp_ele = new Element(temp_p, j);
					temp_n->elements.push_back(temp_ele);
					temp_ng->nodes.push_back(temp_n);
					TotalnumNode++;
					
					
				}
			}			
			else if(numPoly_in_set == 2) { 
				c2++;
				for (int j = 1; j < 5; j++) {
					for (int k = 1; k < 5; k++) {
						Node* temp_n = new Node(this->node_num++);
						//cout << temp_n->idx_adjList << endl;
						poly* temp_p1 = cutline_set[i][0];
						poly* temp_p2 = cutline_set[i][1];
						Element* temp_ele1 = new Element(temp_p1, j);
						Element* temp_ele2 = new Element(temp_p2, k);
						temp_n->elements.push_back(temp_ele1);
						temp_n->elements.push_back(temp_ele2);
						temp_ng->nodes.push_back(temp_n);
						TotalnumNode++;
						//this->node_num++;
						
					}
				}
			}
			else if (numPoly_in_set == 3) {
				c3++;
				for (int j = 1; j < 5; j++) {
					for (int k = 1; k < 5; k++) {
						for (int l = 1; l < 5; l++) {
							Node* temp_n = new Node(this->node_num++);
							//cout << temp_n->idx_adjList << endl;
							poly* temp_p1 = cutline_set[i][0];
							poly* temp_p2 = cutline_set[i][1];
							poly* temp_p3 = cutline_set[i][2];
							Element* temp_ele1 = new Element(temp_p1, j);
							Element* temp_ele2 = new Element(temp_p2, k);
							Element* temp_ele3 = new Element(temp_p3, l);
							temp_n->elements.push_back(temp_ele1);
							temp_n->elements.push_back(temp_ele2);
							temp_n->elements.push_back(temp_ele3);
							temp_ng->nodes.push_back(temp_n);
							TotalnumNode++;
							//this->node_num++;
						}
					}
				}
			}
			else if (numPoly_in_set == 4) {
				c4++;
				for (int j = 1; j < 5; j++) {
					for (int k = 1; k < 5; k++) {
						for (int l = 1; l < 5; l++) {
							for (int m = 1; m < 5; m++) {
								Node* temp_n = new Node(this->node_num++);
								cout << temp_n->idx_adjList << endl;
								poly* temp_p1 = cutline_set[i][0];
								poly* temp_p2 = cutline_set[i][1];
								poly* temp_p3 = cutline_set[i][2];
								poly* temp_p4 = cutline_set[i][3];
								Element* temp_ele1 = new Element(temp_p1, j);
								Element* temp_ele2 = new Element(temp_p2, k);
								Element* temp_ele3 = new Element(temp_p3, l);
								Element* temp_ele4 = new Element(temp_p4, m);
								temp_n->elements.push_back(temp_ele1);
								temp_n->elements.push_back(temp_ele2);
								temp_n->elements.push_back(temp_ele3);
								temp_n->elements.push_back(temp_ele4);
								temp_ng->nodes.push_back(temp_n);
								TotalnumNode++;
								//this->node_num++;
							}
						}
					}
				}
			}
			else {
				cout << "Error: numPoly_in_set > 4" << endl;
			}
			this->nodegroups.push_back(temp_ng);//should be the same length as cutline_set
		}

	}

	//overload to copy vector of polygons
	vector<poly*> operator = (vector<poly*> &vec) {
		vector<poly*> vec_new;
		for (size_t i = 0; i < vec.size(); i++) {
			vec_new.push_back(vec[i]);
		}
		return vec_new;
	}

	row(string row_name, double row_x, double row_y) : name(row_name), x(row_x), y(row_y) {
		this->node_num = 0;
	}

};
//



//for graph
class Graph {
public:
	list<adjNode*> *adjList;
	int *dist;
	int *prev;
	Node* source;
	Node* target;
	row* rowptr;
	int node_Num;

	void add_edge(Node* src, Node* des){// (if all legal)set dest and cal cost => add to list of ajdNode
		//cout << "start to add edge" << endl;

		//check in-node conflict for src
		if (src->elements.size() > 1) {
			for (size_t i = 0; i < src->elements.size(); i++) {
				poly* temp_p1 = src->elements[i]->polygon;
				for (size_t j = i; j < src->elements.size(); j++) {
					if (i == j) continue;
					poly* temp_p2 = src->elements[j]->polygon;

					for (size_t k = 0; k < temp_p1->CG.size(); k++) { //search CG
						if (temp_p1->CG[k] == temp_p2) { //if has in-node conflict
							if (src->elements[i]->color == src->elements[j]->color) {
								//cout << "returned1 in-node conflict for src" << endl;
								return;
							} //and same color...conflict! no edge created
						}
					}
				}
			}
		}


		//check in-node conflict for des
		if (des->elements.size() > 1) {
			for (size_t i = 0; i < des->elements.size(); i++) {
				poly* temp_p1 = des->elements[i]->polygon;
				for (size_t j = i; j < des->elements.size(); j++) {
					if (i == j) continue;
					poly* temp_p2 = des->elements[j]->polygon;

					for (size_t k = 0; k < temp_p1->CG.size(); k++) { //search CG
						if (temp_p1->CG[k] == temp_p2) { //if has in-node conflict

							if (des->elements[i]->color == des->elements[j]->color) {
								//cout << "returned2 in-node conflict for des" << endl;
								return;
							}//and same color...conflict! no edge created
						}
					}
				}
			}
		}

		for(size_t i = 0 ; i< src->elements.size(); i++){
			poly* p_src = src->elements[i]->polygon;
			//check if compatible of elements between two nodes
			for (size_t j = 0; j < des->elements.size(); j++) {
				poly* p_des = des->elements[j]->polygon;

				//same polygon in src, des and not same color -> return
				if (p_src == p_des) {
					if (src->elements[i]->color != des->elements[j]->color) {
						//cout << "returned3" << endl;
						return; }//no edge created
				}

				//diff polygons in src, des and same color and in CG =>return
				else {
					for (size_t k = 0; k < p_src->CG.size(); k++) { //search CG
						if (p_src->CG[k] == p_des) { //if has conflict
							if (src->elements[i]->color == des->elements[j]->color) { 
								//cout << "returned4" << endl;
								return; } //and same color...conflict! no edge created
						}
					}
				}
			}

		}

		int E_num = 0;
		for (size_t i = 0; i < des->elements.size(); i++) {
			if (des->elements[i]->color == 4) E_num++;
		}
		
		//cout << "added" << endl;
		TotalnumEdge++;
		adjNode* new_aNode = new adjNode(src, des, E_num);
		//cout <<"src->idx_adjList: "<< src->idx_adjList << endl;
		adjList[src->idx_adjList].push_back(new_aNode);
	}

	void add_edge_for_start_node() {
		NodeGroup* temp_ng = this->rowptr->nodegroups[0]; //source link to the first set 
		for (size_t i = 0 ; i < temp_ng->nodes.size(); i++) {
			int E_num = 0;
			Node* temp_des = temp_ng->nodes[i];
			for (size_t j = 0; j < temp_des->elements.size(); j++) {
				if (temp_des->elements[j]->color == 4) E_num++;
			}
			adjNode* new_aNode = new adjNode(source, temp_des, E_num);
			adjList[source->idx_adjList].push_back(new_aNode);
			TotalnumEdge++;
		}
		

	}

	void add_edge_for_sink_node() {
		int size_ng = this->rowptr->nodegroups.size();
		NodeGroup* temp_last_ng = this->rowptr->nodegroups[size_ng-1];
		for (size_t i = 0; i < temp_last_ng->nodes.size(); i++) {
			Node* temp_src = temp_last_ng->nodes[i];
			adjNode* new_aNode = new adjNode(temp_src, target, 0);
			adjList[target->idx_adjList].push_back(new_aNode);
			TotalnumEdge++;
		}

	}

	int get_smallest_idx_inQ(list<int> Q) {  //...time comp: O(N)
		list<int> ::iterator it;
		int smallest_dist = infinity;
		int smallest_idx = -1;
		for (it = Q.begin(); it != Q.end(); it++) {
			if (this->dist[*it] < smallest_dist) { 
				smallest_dist = this->dist[*it];
				smallest_idx = *it;
			}
		}
		return smallest_idx;
	}

	//int get_cost_inadjList(list<adjNode*> *aList, int picked, int dest) {
	//	list<adjNode*> temp_adList_des = aList[dest];
	//	list<adjNode*> ::iterator it;
	//	for (it = aList.begin(); it != aList.end(); it++) {
	//		if (this->dist[*it] < smallest_dist) {
	//			smallest_dist = this->dist[*it];
	//			smallest_idx = *it;
	//		}
	//	}
	//	return smallest_idx;
	//}


	void solve_dijk() {
		int start_idx = source->idx_adjList;

		list<int> Q;
		set<int> S;

		this->dist[start_idx] = 0;

		for (int i = 0; i < this->node_Num + 2; i++) {//skip placing source node in Q
			if(!(i == start_idx)) Q.push_back(i);
			//cout << this->dist[i] << endl;
		}
		S.insert(start_idx);
		list<adjNode*> ::iterator it_;

		

		//do the start node
		for (it_ = this->adjList[start_idx].begin(); it_ != this->adjList[start_idx].end(); it_++) {
			adjNode* temp_aNode = *it_;
			Node* dest_Node = temp_aNode->next;
			int dest_idx_adjList = dest_Node->idx_adjList;
			//cout << (this->dist[start_idx] + (temp_aNode->cost)) << ' ' << this->dist[dest_idx_adjList] << endl;
			//cout << start_idx << "  " << dest_idx_adjList << endl;
			if ((this->dist[start_idx] + (temp_aNode->cost)) < this->dist[dest_idx_adjList]) { //relax (u,v)
				this->dist[dest_idx_adjList] = (this->dist[start_idx] + (temp_aNode->cost));
				this->prev[dest_idx_adjList] = start_idx;
				//cout << "updated dist: "<< this->dist[dest_idx_adjList] << endl;
			}
		}
		
		//do the rest nodes
		while (!Q.empty()) {
			int idx_src = this->get_smallest_idx_inQ(Q);
			//cout <<"idx_src" <<idx_src << endl;
			list<adjNode*> ::iterator it;
			if (idx_src == -1) break;
			for (it = this->adjList[idx_src].begin(); it != this->adjList[idx_src].end(); it++) {
				//cout << "hello0" << endl;
				adjNode* temp_aNode = *it;
				//cout << "hello1" << endl;
				Node* dest_Node = temp_aNode->next;
				//cout << "hello2" << endl;
				int dest_idx_adjList = dest_Node->idx_adjList;
				//cout << "hello3" << endl;
				//cout << (this->dist[start_idx] + (temp_aNode->cost)) << ' ' << this->dist[dest_idx_adjList] << endl;
				//cout << start_idx << "  " << dest_idx_adjList << endl;
				if ((this->dist[idx_src] + (temp_aNode->cost)) < this->dist[dest_idx_adjList]) { //relax (u,v)
					this->dist[dest_idx_adjList] = (this->dist[idx_src] + (temp_aNode->cost));
					this->prev[dest_idx_adjList] = idx_src;
					//cout << "updated dist: "<< this->dist[dest_idx_adjList] << endl;
				}
			}
			Q.remove(idx_src);
			//cout << "removed from Q " << idx_src << endl;
			S.insert(idx_src);
		}


	}

	int find_last_node_idx() {
		row* temp_row = this->rowptr;

		for(unsigned int i = temp_row->nodegroups.size()-1; i > 0 ; i--){
			for (unsigned int j = 0; j < temp_row->nodegroups[i]->nodes.size(); i++) {
				int idx_node = temp_row->nodegroups[i]->nodes[j]->idx_adjList;
				if (this->dist[idx_node] < 100000000) { 
					//cout << this->dist[idx_node]<<" idx: "<< idx_node << endl;

					//cout << "idx_adjList: " << temp_aNode->key->idx_adjList << endl;
					return idx_node;
				}
			}
		}
		return -1;

		//int idx_target = this->target->idx_adjList;
		//int idx_source = this->source->idx_adjList;

		//for (int i = this->node_Num; i != 0; i--) {
		//	if (this->dist[i] < 100000000) {
		//		cout << this->dist[i] << endl;
		//	}

		//}
		////bool flag = false;
		////if ((this->dist[idx_target]) < 100000000) flag = true;
		//cout << (this->dist[idx_target]) <<" "<< this->prev[idx_target] << endl;
		////cout << (this->dist[idx_source+1]) <<" "<< this->prev[idx_source+1] << endl;
	
	}

	void assign_colors() {
		int cur_idx = this->find_last_node_idx();
		if (cur_idx == -1) return;
		while (cur_idx != source->idx_adjList) {
			int prev_idx = this->prev[cur_idx];

			adjNode* temp_aNode = *(this->adjList[cur_idx].begin());
			Node* temp_n = temp_aNode->key;
			for (unsigned int i = 0; i < temp_n->elements.size(); i++) {
				poly* temp_p = temp_n->elements[i]->polygon;
				int color_idx = temp_n->elements[i]->color;
				char color;
				if (color_idx == 1) color = 'A';
				else if (color_idx == 2) color = 'B';
				else if (color_idx == 3) color = 'C';
				else color = 'E';

				temp_p->color = color;
			}
			cur_idx = prev_idx;
		}


	}

	//int n;
	Graph(int nodeNum, row* myrow) {
		adjList = new list<adjNode*>[nodeNum + 2]; //two for start/end node
		//nodeList = new list<Node*>[nodeNum + 2]; //two for start/end node
		dist = new int[nodeNum + 2];
		prev = new int[nodeNum + 2];
		for (int i = 0; i < nodeNum + 2; i++) {
			dist[i] = infinity;
			prev[i] = -1;
		}

		//adjNode* source_temp = new adjNode(nodeNum); //Place the start node to the last position-1
		Node* source_temp = new Node(nodeNum); //Place the start node to the last position-1
		source = source_temp;
		//
		//adjNode* target_temp = new adjNode(nodeNum+1); //Place the end node to the last position
		Node* target_temp = new Node(nodeNum+1); //Place the end node to the last position
		target = target_temp;

		this->rowptr = myrow;
		this->node_Num = nodeNum;
		//target_temp->next = nullptr;
		//this->n = nodeCount;
		//this->adjList = new list<Node*>[n];
	}
};
//

//global var 2
vector<row*> rows;
vector<Graph*> graphs;
//

void print_personal_info() {
	cout << "===========================================================" << endl;
	cout << "                   VLSI DFM Final Project                  " << endl;
	cout << "                    Date: 23th June 2021                   " << endl;
	cout << "                    Name: Weishiun HUNG                    " << endl;
	cout << "       School: National Tsing Hua University, Taiwan       " << endl;
	cout << "                   Student ID: 109065527                   " << endl;
	cout << "===========================================================" << endl;
}

void parcer(char* filename) {
	cout << "input file is: " << filename << endl;

	ifstream infile(filename);
	string line;

	//ROW

	getline(infile, line);
	istringstream iss_row_start(line);
	string linehead;
	iss_row_start >> linehead;
	while (linehead == "ROW")
	{
		istringstream iss_row(line);
		string row_name;
		double row_x, row_y;
		iss_row >> linehead >> row_name >> row_x >> row_y;
		row* R = new row(row_name, row_x, row_y);
		//R->name = row_name;
		rows.push_back(R);
		//cout << linehead <<"\t"<<row_name<<"\t"<< row_x << "\t" << row_y << endl;
		//POLYGON
		getline(infile, line);
		istringstream iss_poly_start(line);
		iss_poly_start >> linehead;
		while (linehead == "POLYGON") {

			istringstream iss_poly(line);
			string poly_name;
			num_poly++; // count number of polygons
			iss_poly >> linehead >> poly_name;

			//cout << "\t"<<linehead << "\t" << poly_name  << endl;
			//MASK
			getline(infile, line);
			istringstream iss_mask_start(line);
			char mask_name;
			iss_mask_start >> linehead >> mask_name;
			poly* P = new poly(mask_name, poly_name);
			//P(mask_name, poly_name);
			R->polygons.push_back(P);
			//cout << "\t" << "\t" << linehead <<"\t"<< mask_name << endl;			
			//LAYER
			getline(infile, line);
			//cout << line << endl;
			//RECT
			getline(infile, line);
			istringstream iss_rect_start(line);
			iss_rect_start >> linehead;
			while (linehead == "RECT") {
				istringstream iss_rect(line);
				double x_left, y_left, x_right, y_right;
				iss_rect >> linehead >> x_left >> y_left >> x_right >> y_right;
				rect* RE = new rect(x_left, y_left, x_right, y_right);
				P->push_rect(RE);
				//cout << "\t" << "\t" << "\t" << linehead  << "\t" << x_left<< "\t" << y_left <<"\t"<< x_right<<"\t"<< y_right << endl;
				//next line
				getline(infile, line); // (RECT X X X X )or (END)
				istringstream iss_rect_Linehead_detect(line);
				iss_rect_Linehead_detect >> linehead;
			}
			//cout << line << endl; //END
			//
			getline(infile, line); //END P0
			//cout << line << endl;
			//
			getline(infile, line);
			istringstream iss_poly_Linehead_detect(line);
			iss_poly_Linehead_detect >> linehead;
		}
		//cout << line << endl;
		getline(infile, line);
		istringstream iss_row_Linehead_detect(line);
		iss_row_Linehead_detect >> linehead;

		// process pair (a,b)
	}
}

void outputfile(char* filename) {
	cout << "output file is: " << filename << endl;
	ofstream outfile(filename);
	for (size_t i = 0; i < rows.size(); i++) {
		outfile << "ROW " << rows[i]->name << " " << rows[i]->x << " " << rows[i]->y << " " << endl;
		for (size_t j = 0; j < rows[i]->polygons.size(); j++) {
			poly* temp_p = rows[i]->polygons[j];
			outfile << "\t" << "POLYGON " << temp_p->name << endl;
			outfile << "\t" << "\t" << "MASK " << temp_p->color << endl;
			outfile << "\t" << "\t" << "LAYER M1" << endl;
			for (size_t k = 0; k < rows[i]->polygons[j]->rects.size(); k++) {
				rect* temp_rect = rows[i]->polygons[j]->rects[k];
				outfile << "\t" << "\t" << "\t" << "RECT " << temp_rect->xl << " " << temp_rect->yl << " " << temp_rect->xr << " " << temp_rect->yr << endl;

			}
			outfile << "\t" << "\t" << "END" << endl;
			outfile << "\t" << "END " << temp_p->name << endl;
		}
		outfile << "END " << rows[i]->name << endl;
	}

}


//test functions
void print_left_boundaries(row* R) {
	for (unsigned int i = 0; i < R->polygons.size(); i++) {
		cout << R->polygons[i]->min_x << " ";
	}
	cout << endl;
}
//

void build_graph() {}
void solve_dijkstra() {}
void assign_color() {}

int main(int argc, char* argv[]) {
	auto t1 = chrono::system_clock::now();
	print_personal_info();
	dmin = atof(argv[1]);
	cout << "dmin specified= " << dmin << endl;
	parcer(argv[2]);
	cout << "#polygons: " << num_poly << endl;
	auto t2 = chrono::system_clock::now();
	//
	//print_left_boundaries(rows[0]);

	//initialization
	auto t_c1 = chrono::system_clock::now();

	for (unsigned int i = 0; i < rows.size(); i++) {
		rows[i]->build_cg();
		//cout << "h4" << endl;
		rows[i]->set_cutline();
		//cout << "h5" << endl;
		rows[i]->set_dummy_ext_row();
		rows[i]->create_cutline_set();
		//cout << rows[i]->polygons.size() << " " << rows[i]->sorted_poly.size() << " " << rows[i]->cutline_def.size() << " " << rows[i]->cutline_set.size() << endl;
		
	}
	
	auto t_c2 = chrono::system_clock::now();
	auto t_c3 = chrono::system_clock::now();

	if (numAs == num_poly) {//skip graph and output result
		cout << "TotalnumNode: " << TotalnumNode << endl;
		cout << "#iso-polygons: " << numAs << " (all painted A)" << endl;
	}
	else {

		for (unsigned int i = 0; i < rows.size(); i++) {
			rows[i]->create_nodes();

			Graph* new_graph = new Graph(rows[i]->node_num, rows[i]);

			//start point to first set of nodes
			new_graph->add_edge_for_start_node();

			//create rest of the edges
			for (size_t j = 0; j < rows[i]->nodegroups.size() - 1; j++) {
				NodeGroup* temp_ng_src = rows[i]->nodegroups[j];
				NodeGroup* temp_ng_des = rows[i]->nodegroups[j + 1];

				for (size_t k = 0; k < temp_ng_src->nodes.size(); k++) {
					for (size_t l = 0; l < temp_ng_des->nodes.size(); l++) {
						new_graph->add_edge(temp_ng_src->nodes[k], temp_ng_des->nodes[l]);

					}
				}
			}
			//last set to sink node
			new_graph->add_edge_for_sink_node();
			graphs.push_back(new_graph); //complete graph

			//ccccc = ccccc + rows[i]->sorted_poly.size() - rows[i]->cutline_def.size();
			//cout << ccccc << endl;

		}
		

		t_c2 = chrono::system_clock::now();
		t_c3 = chrono::system_clock::now();
		//solve graph start
		for (unsigned int i = 0; i < graphs.size(); i++) {
			graphs[i]->solve_dijk();
			graphs[i]->assign_colors();
		}
		//cout << "c1, c2, c3, c4: " << c1 << ' ' << c2 << ' ' << c3 << ' ' << c4 << endl;
		//cout << "TotalnumNode: " << TotalnumNode << endl;
		//cout << "TotalnumEdge: " << TotalnumEdge << endl;
		//cout << "#iso-polygons: " << numAs << " (all painted A)" << endl;


	}

	//
	auto t_c4 = chrono::system_clock::now();
	//cal ends
	auto t3 = chrono::system_clock::now();
	outputfile(argv[3]);
	auto t4 = chrono::system_clock::now();

	cout << "IO time: " << chrono::duration_cast<chrono::duration<double>>(t2 - t1 + t4 - t3).count() << endl;
	cout << "Initialization time: " << chrono::duration_cast<chrono::duration<double>>(t_c2 - t_c1).count() << endl;
	cout << "Solving time: " << chrono::duration_cast<chrono::duration<double>>(t_c4 - t_c3).count() << endl;
	cout <<"Total runtime: " << chrono::duration_cast<chrono::duration<double>>(t4- t1).count() << endl;

	return 0;
}