#include "all_headers.hpp"
#include "TaxTree.hpp"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

#undef DEBUGME
//#define DEBUGME

using namespace std;
using namespace metag;


template <class tid_T>
void TaxTree<tid_T>::printPathToRoot(tid_T tid, bool with_names) {
  vector<tid_T> t;
  getPathToRoot(tid, t);
  if (with_names) {
    cout << tid  << " ";
    typename TaxTree<tid_T>::const_iterator t = find(tid);
    if (t == (*this).end()) {
      cout << " TaxNode not found! ";
    } else {
      cout << t->second->name() << " ";
    }
    cout << " :: ";
  } else {
    cout << tid << " :: ";
  }
  for (size_t j=0; j<t.size(); j++) {
    if (with_names) {
      cout << t[j] << " ";
      typename TaxTree::const_iterator t = find(tid);
      if (t == (*this).end()) {
        cout << " TaxNode not found! ";
      } else {
        cout << t->second->name() << " ";
      }
    } else {
      cout << t[j] << " ";
    }
  }
  cout << endl;
}

template<class tid_T>
void TaxTree<tid_T>::reversePathsToRoot() {
  for (typename TaxTree::iterator t = (*this).begin(); t != (*this).end(); t++) {
    t->second->reversePathToRoot();
  }
}


/*
template<class tid_T>
void TaxTree<tid_T>::getPathToRoot(tid_T tid, vector<tid_T> &out) {
  out.clear();

  //sanity check: a TaxNode for tid exists
  if ((*this).find(tid) == (*this).end()) {
    registerFailure(tid);
    return;
  }

  TaxNode<tid_T> *current = (*this).find(tid)->second;
  if (current->parent() == current->id()) {
    return;
  }

  if ((*this).find(current->parent()) == (*this).end()) {
    cerr << "failed to find parent TaxNode for taxid " << current->id() << " whose parent is " << current->parent() << endl;
    cerr << "fatal error!\n";
    exit(-1);
  }
  current = (*this).find(current->parent())->second;
  out.push_back(current->id());

  while (true) {
    //test for root
    if (current->parent() == current->id()) {
      return;
    }

    current = (*this).find(current->parent())->second;
    out.push_back(current->id());
  }
}
*/

//returns LCA; on return, children will contain all entries in the input tax_ids set, plus the lca
//plus all tax IDs between the tax_ids and the LCA
/*
template<class tid_T>
tid_T TaxTree<tid_T>::getLcaMap(const set<tid_T> &tax_ids, __gnu_cxx::hash_map<tid_T, set<tid_T> > &children) const {
  children.clear();

 // children[1]; //add root
  set<tid_T> good_tax_ids;

  //add to the set all tax IDs between the leaves (tax_ids) and root
  for (typename set<tid_T>::iterator t = tax_ids.begin(); t != tax_ids.end(); t++) {
    if (this->find(*t) != this->end()) {
      TaxNode<tid_t> *tn = this->find(*t)->second;
      if (!tn) {
        registerErasedTid(*t);
        continue;
      }  
      good_tax_ids.insert(*t);
      assert(tn); //should be impossible to trip ...

      //get path to root; not that the path does not include *t
      const vector<tid_T> *p = tn->getPathToRoot();
      assert(p);

      if (p->size()) {
        children[(*p)[0]].insert(*t);
      }  
      for (size_t j=1; j<p->size(); j++) {
        assert(this->find((*p)[j]) != this->end());
        children[(*p)[j]].insert((*p)[j-1]);
      }
    }
  }

  assert(children.find(1) != children.end());

  //delete nodes between root and the LCA; a node can be deleted
  //if it has a single child, wrt the children map
  tid_T lca = 1;
  vector<tid_T> eraseme;
  while (true) {
    if (children[lca].size() == 1) {
      tid_T child = *(children[lca].begin());
      if (tax_ids.find(child) == tax_ids.end()) {
        eraseme.push_back(lca);
        lca = child;
      } else {
        break;
      }
    } else {
      break;
    }
  }
  for (size_t j=0; j<eraseme.size(); j++) {
    children.erase(eraseme[j]);
  }

  //add in the starting tid set
  for (typename set<tid_T>::const_iterator t = good_tax_ids.begin(); t != good_tax_ids.end(); t++) {
    children[*t];
  }

  return lca;
}
*/


/*
template<class tid_T>
void TaxTree<tid_T>::setRanks(const char *fn) {
  ifstream in(fn);
  assert(in);
  map<tid_T, string> mp;
  tid_T tid;
  string rank;
  while (in >> tid >> rank) {
    mp[tid] = rank;
  }

  string parent_rank;
  for (typename TaxTree::iterator t = this->begin(); t != this->end(); t++) {
    rank = ""; 
    parent_rank = "";
    tid_T id = t->second->id() ;
    tid_T parent = t->second->parent();

    if (mp.find(parent) != mp.end()) {
      parent_rank = mp[parent];
    } else {
      cout << "failed to find rank for parent taxid: " << parent << endl;
    }

    if (mp.find(id) != mp.end()) {
      rank = mp[id];
    } else {
      cout << "failed to find rank for taxid: " << id << endl;
    }
    t->second->setRank(rank, parent_rank);
  }
}
*/

