#ifndef __TAX_TREE__
#define __TAX_TREE__

#include <ext/hash_map>
#include <set>
#include <map>
#include <iostream>
#include "TaxNode.hpp"

namespace metag {

/**
*/

//typedef __gnu_cxx::hash_map<uint32_t, TaxNode*> TaxNodeHash;


template<class tid_T>
class TaxTree : public __gnu_cxx::hash_map<tid_T, TaxNode<tid_T>*> {
public :

  TaxTree(const char *nodes_fn) {

    ifstream in(nodes_fn);
    if (! in.is_open()) {
      cerr << "failed to open " << nodes_fn << " for reading\n";
      exit(-1);
    }

    string line;
    getline(in, line);
    getline(in, line);
    int count;
    //the count may not be accurate, so ignore it
    in >> count;
    getline(in, line);
    TaxNode<tid_T> *t_last = 0;
    int ct = 0;
    while (true) {
      ++ct;
      streampos p = in.tellg();
      if (in.eof()|| ! in.good() || (int)p == -1) {
        break;
      }
      TaxNode<tid_T> *t = TaxNode<tid_T>::read(in);
      t_last = t;
      (*this)[t->id()] = t;
    }
    in.close();

    for (typename TaxTree::iterator t = (*this).begin(); t != (*this).end(); t++) {
      t->second->setPathToRoot(this);
    }
    (*this)[1]->eraseChild(1);
  }

  //void getPathToRoot(tid_T tid, std::vector<tid_T> &out);
  void getPathToRoot(tid_T tid, vector<tid_T> &out) {
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

  string getName(tid_T id) const {
    return ( this->find(id) == this->end() ) ? "" : this->find(id)->second->name();
  } 
  void printChildren(tid_T id) const {
    const std::set<tid_T> &s = this->find(id)->second->getChildren();
    std::cout << id << " :: ";
    for (typename std::set<tid_T>::const_iterator t3 = s.begin(); t3 != s.end(); t3++) {
      std::cout << *t3 << " ";
    }
    std::cout << std::endl;
  }

  void reversePathsToRoot();

  void setRanks(const char *fn) {
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

  void printPathToRoot(tid_T tid, bool with_names);

  void registerErasedTid(tid_T tid) const {
    if (m_erased.find(tid) == m_erased.end()) {
      std::cerr << "from TaxTree.hpp: erased tid " << tid << " from tax ID input set, since tid not found in tree; nodes in tree: " << this->size() << "\n";
      m_erased.insert(tid);
    }
  }

  void registerFailure(tid_T tid, const char *msg = 0) const {
    if (m_failed.find(tid) == m_failed.end()) {
      std::cerr << "from TaxTree.hpp: failed to find taxid " << tid << " in tax_tree; nodes in tree: " << this->size() << "\n";
      if (msg) {
        cout << " :: " << msg << std::endl;
      }
      m_failed.insert(tid);
    }
  }

  //! returns LCA; on return, children will contain all entries in the input tax_ids set, plus the lca
  //! plus all tax IDs between the tax_ids and the LCA
  tid_T getLcaMap(const set<tid_T> &tax_ids, __gnu_cxx::hash_map<tid_T, set<tid_T> > &children) const {
    children.clear();

    if (! tax_ids.size()) {
      cerr << "from TaxTree.hpp: error: input tax_ids set is empty\n";
      exit(-1);
    }

    if (tax_ids.size() == 1) {
      typename set<tid_T>::const_iterator t = tax_ids.begin();
      tid_T tid = *t;
      if (this->find(tid) == this->end()) {
        registerErasedTid(tid);
        tid_T x = ~0;
        return x;
      }
      children[tid];
      return tid;
    }

    set<tid_T> good_tax_ids;

    //add to the 'children' set all tax IDs between the leaves (tax_ids) and root
    for (typename set<tid_T>::iterator t = tax_ids.begin(); t != tax_ids.end(); t++) {
      tid_T cur_tid = *t;
      if (this->find(cur_tid) != this->end()) {
        TaxNode<tid_T> *tn = this->find(cur_tid)->second;
        if (!tn) {
          registerErasedTid(cur_tid);
          continue;
        }
        good_tax_ids.insert(cur_tid);
        assert(tn); //n.b. should be impossible to trip ...

        //get path to root for the cur_tid;  not that the path does not include cur_tid
        const vector<tid_T> *p = tn->getPathToRoot();
        assert(p);
        assert(p->back() == 1);

        //add cur_tid as a child of its parent
        if (p->size()) {
          children[(*p)[0]].insert(cur_tid);
        }

        assert(p->size() > 1);
        for (size_t j=1; j<p->size(); j++) {
          assert(this->find((*p)[j]) != this->end());
          children[(*p)[j]].insert((*p)[j-1]);
        }
      } else {
         registerErasedTid(cur_tid);
      }
    }

    if (!good_tax_ids.size()) {
      cerr << "\nfrom TaxTree.hpp: WARNING: good_tax_ids.size() is empty\n";
      for (typename set<tid_T>::const_iterator t2 = tax_ids.begin(); t2 != tax_ids.end(); t2++) {
        if (this->find(*t2) == this->end()) {
          cerr << "    input tax id " << *t2 << " was NOT found in the taxtree\n";
        } else {
          cerr << "    input tax id " << *t2 << " WAS found in the taxtree (something's wrong if you're seeing this!)\n";
        }
      }
      cerr << endl;
      return 0;
    }

    if (children.find(1) == children.end()) {
      cerr << "from TaxTree.hpp: failed to find root in 'children' map.\n";
      cerr << "map.size(): " << children.size() << endl;
      exit(-1);
    }

    //delete nodes between root and the LCA; a node can be deleted
    //if it has a single child, wrt the children map
    tid_T lca = 1;
    vector<tid_T> eraseme;
    while (true) {
      if (children[lca].size() == 1) {
        tid_T child = *(children[lca].begin());
        if (tax_ids.find(lca) == tax_ids.end()) {
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

private:

  std::map<tid_T, tid_T> m_lca_options;

  mutable std::set<tid_T> m_failed;
  mutable std::set<tid_T> m_erased;

};

}

#endif
