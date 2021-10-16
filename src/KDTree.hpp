// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_
#include <exception>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <queue>
#include "Point.hpp"
using namespace std;
//nodo
template <size_t N>
struct KDTreeNode {

    Point<N> val;
    size_t arrpos;
    KDTreeNode<N>* nodos[2] = { nullptr,nullptr };

    KDTreeNode(const Point<N>& value,size_t arrpos);

};
template<size_t N>
KDTreeNode<N>::KDTreeNode(const Point<N>& value, size_t arrpos)
{
    val = value;
    this->arrpos = arrpos;
}




template <size_t N, typename ElemType>
class KDTree {
public:
    typedef std::pair<Point<N>, ElemType> value_type;

    KDTree();

    ~KDTree();

    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);

    size_t dimension() const;

    size_t size() const;
    bool empty() const;

    bool contains(const Point<N>& pt) const;
    bool find(const Point<N>& pt,KDTreeNode<N>**& p);
    void insert(const Point<N>& pt, const ElemType& value);

    ElemType& operator[](const Point<N>& pt);

    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;

    ElemType knn_value(const Point<N>& key, size_t k) const;

    std::vector<ElemType> knn_query(const Point<N>& key, size_t k) const;

private:



    KDTreeNode<N>* root;
    ElemType d;
    size_t dimension_;
    size_t size_;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    dimension_ = N;
    size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
    size_ = rhs.size_;
    dimension_ = rhs.dimension_;
    root = rhs.root;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {
    KDTree<N, ElemType>* t=new  KDTree<N, ElemType>();
    queue<KDTreeNode<N>*>cola;
    cola.push(rhs.root);
    while (!cola.empty()) {
        KDTreeNode<N>* top=cola.front();
        if (top->nodos[0]) cola.push(top->nodos[0]);
        if (top->nodos[0]) cola.push(top->nodos[1]);
        (*t).insert(top->val,top->arrpos);
        cola.pop();
    }
    t->size_ = rhs.size_;
    t->dimension_ = dimension_;
    return *t;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    // TODO(me): Fill this in.
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    
    return (size_==0);
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
    if (root == nullptr) {
        return false;
    }
    else {
        KDTreeNode<N>** p= const_cast<KDTreeNode<N>**>(&root);

        while (*p && (*p)->val != pt) {
            bool opt = pt > (*p)->val;
            p = &((*p)->nodos[opt]);
        }
        return (*p != 0 && (*p)->val == pt);
    }
}
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(const Point<N>& pt, KDTreeNode<N>**& p) {
    if (root == nullptr) {
        return false;
    }
    else {
        while (*p && (*p)->val != pt) {
            bool opt = pt > (*p)->val;
            p = &((*p)->nodos[opt]);
        }
       
        return (*p != 0 && (*p)->val==pt);
    }
}
template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
    if (root == nullptr) {
        root= new KDTreeNode<N>(pt,value);
    }
    else {
        KDTreeNode<N>** p = const_cast<KDTreeNode<N>**>(&root);
        if (find(pt, p)) {
            (*p)->arrpos=value;
            size_--;
        }
        else {
            *p = new KDTreeNode<N>(pt, value);
        }
        
    }
    size_++;
    
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
    KDTreeNode<N>** p = const_cast<KDTreeNode<N>**>(&root);
    if (find(pt, p)) {
        return (*p)->arrpos;
    }
    else {
        *p = new KDTreeNode<N>(pt, 0);
        size_++;
        return (*p)->arrpos;
    }
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt){
    if (root == nullptr) {
        throw  out_of_range("out_of_range");
    }
    else {
        KDTreeNode<N>** p = const_cast<KDTreeNode<N>**>(&root);
        while (*p && (*p)->val != pt) {
            bool opt = pt > (*p)->val;
            p = &((*p)->nodos[opt]);
        }

        if (p != nullptr) {
            return (*p)->arrpos;
        }
        else {
            throw  out_of_range("out_of_range");
        }
    }
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {

    if (root == nullptr) {
        throw  out_of_range("out_of_range");
    }
    else {
        KDTreeNode<N>** p = const_cast<KDTreeNode<N>**>(&root);
        while (*p && (*p)->val != pt) {
            bool opt = pt > (*p)->val;
            p = &((*p)->nodos[opt]);
        }

        if (p != nullptr) {
            return (*p)->arrpos;
        }
        else {
            throw  out_of_range("out_of_range");
        }
    }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N>& key, size_t k) const {
    // TODO(me): Fill this in.
    ElemType new_element;
    return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N>& key,size_t k) const {
    // TODO(me): Fill this in.
    std::vector<ElemType> values;
    return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class


#endif  // SRC_KDTREE_HPP_


