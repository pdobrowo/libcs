/**
 * Copyright (C) 2009-2015  Przemysław Dobrowolski
 *
 * This file is part of the Configuration Space Library (libcs), a library
 * for creating configuration spaces of various motion planning problems.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef LIBCS_PREDICATE_SENTENCE_3_H
#define LIBCS_PREDICATE_SENTENCE_3_H

#include "Spin_3.h"
#include <CGAL/Kernel/global_functions.h>
#include <memory>
#include <cassert>

namespace CS
{
template<class Kernel_>
class Spin_quadric_tree_3;

template<class Kernel_>
Spin_quadric_tree_3<Kernel_> operator &(const Spin_quadric_tree_3<Kernel_> &left, const Spin_quadric_tree_3<Kernel_> &right);

template<class Kernel_>
Spin_quadric_tree_3<Kernel_> operator |(const Spin_quadric_tree_3<Kernel_> &left, const Spin_quadric_tree_3<Kernel_> &right);

template<class Kernel_>
class Spin_quadric_tree_3
{
    typedef typename Kernel_::RT             RT;
    typedef typename Kernel_::Spin_3         Spin_3;
    typedef typename Kernel_::Vector_3       Vector_3;

    typedef typename Kernel_::Spin_quadric_3 Spin_quadric_3;

    class Node;
    typedef boost::shared_ptr<Node> NodePtr;

    class Node
    {
    public:
        struct Operator
        {
            enum class Type
            {
                And,
                Or
            };
        };

        NodePtr clone() const
        {
            NodePtr node(new Node());

            node->m_leaf = m_leaf;

            if (m_leaf)
            {
                node->m_quadric = m_quadric;
            }
            else
            {
                node->m_operator = m_operator;
                node->m_left = m_left->clone();
                node->m_right = m_right->clone();
            }

            return node;
        }

        static NodePtr leaf(const Spin_quadric_3 &quadric)
        {
            NodePtr node(new Node());

            node->m_leaf = true;
            node->m_quadric = quadric;

            return node;
        }

        static NodePtr join(const NodePtr &left, const NodePtr &right, typename Operator::Type op)
        {
            NodePtr node(new Node());

            node->m_leaf = false;
            node->m_operator = op;
            node->m_left = left;
            node->m_right = right;

            return node;
        }

        size_t size() const
        {
            if (m_leaf)
                return 1;

            return m_left->size() + m_right->size();
        }

        bool evaluate(const Spin_3 &spin)
        {
            // FIXME: This is naive; use Canny's O(lg(n)) implementation
            if (m_leaf)
                return m_quadric.evaluate_sign(spin) != CGAL::NEGATIVE;

            switch (m_operator)
            {
            case Operator::And:
                return m_left->evaluate(spin) && m_right->evaluate(spin);

            case Operator::Or:
                return m_left->evaluate(spin) || m_right->evaluate(spin);

            default:
                assert(0);
            }
        }

    private:
        bool                    m_leaf;

        // leaf
        Spin_quadric_3          m_quadric;

        // node
        typename Operator::Type m_operator;
        NodePtr                 m_left;
        NodePtr                 m_right;
    };

    NodePtr         m_root;

    // construct
    Spin_quadric_tree_3(NodePtr root);

public:
    // construct empty sentence
    Spin_quadric_tree_3();

    // construct a sentence of one quadric
    explicit Spin_quadric_tree_3(const Spin_quadric_3 &quadric);

    // copy and construct
    Spin_quadric_tree_3(const Spin_quadric_tree_3 &other);
    Spin_quadric_tree_3 &operator =(const Spin_quadric_tree_3 &other);

    // evaluation
    bool evaluate(const Spin_3 &spin);

    // iteration
    //const_iterator begin() const;
    //const_iterator end() const;

    // container
    size_t size() const;

    // operators
    //Spin_quadric_tree_3 &operator |= (const Spin_quadric_tree_3 &other);
    //Spin_quadric_tree_3 &operator &= (const Spin_quadric_tree_3 &other);

    friend Spin_quadric_tree_3<Kernel_> operator & <>(const Spin_quadric_tree_3<Kernel_> &left, const Spin_quadric_tree_3<Kernel_> &right);
    friend Spin_quadric_tree_3<Kernel_> operator | <>(const Spin_quadric_tree_3<Kernel_> &left, const Spin_quadric_tree_3<Kernel_> &right);
};

template<class Kernel_>
Spin_quadric_tree_3<Kernel_>::Spin_quadric_tree_3(NodePtr root)
    : m_root(root)
{
}

template<class Kernel_>
Spin_quadric_tree_3<Kernel_>::Spin_quadric_tree_3()
{
}

template<class Kernel_>
Spin_quadric_tree_3<Kernel_>::Spin_quadric_tree_3(const Spin_quadric_3 &quadric)
{
    m_root = Node::leaf(quadric);
}

template<class Kernel_>
Spin_quadric_tree_3<Kernel_>::Spin_quadric_tree_3(const Spin_quadric_tree_3 &other)
{
    if (other.m_root)
        m_root = other.m_root->clone();
}

template<class Kernel_>
Spin_quadric_tree_3<Kernel_> &Spin_quadric_tree_3<Kernel_>::operator =(const Spin_quadric_tree_3 &other)
{
    if (other.m_root)
        m_root = other.m_root->clone();

    return *this;
}

template<class Kernel_>
Spin_quadric_tree_3<Kernel_> operator |(const Spin_quadric_tree_3<Kernel_> &left, const Spin_quadric_tree_3<Kernel_> &right)
{
    if (!left.m_root && !right.m_root)
        return Spin_quadric_tree_3<Kernel_>();

    if (!right.m_root)
        return Spin_quadric_tree_3<Kernel_>(left.m_root->clone());

    if (!left.m_root)
        return Spin_quadric_tree_3<Kernel_>(right.m_root->clone());

    return Spin_quadric_tree_3<Kernel_>(
                Spin_quadric_tree_3<Kernel_>::Node::join(
                    left.m_root->clone(),
                    right.m_root->clone(),
                    Spin_quadric_tree_3<Kernel_>::Node::Operator::Or));
}

template<class Kernel_>
Spin_quadric_tree_3<Kernel_> operator &(const Spin_quadric_tree_3<Kernel_> &left, const Spin_quadric_tree_3<Kernel_> &right)
{
    if (!left.m_root && !right.m_root)
        return Spin_quadric_tree_3<Kernel_>();

    if (!right.m_root)
        return Spin_quadric_tree_3<Kernel_>(left.m_root->clone());

    if (!left.m_root)
        return Spin_quadric_tree_3<Kernel_>(right.m_root->clone());

    return Spin_quadric_tree_3<Kernel_>(
                Spin_quadric_tree_3<Kernel_>::Node::join(
                    left.m_root->clone(),
                    right.m_root->clone(),
                    Spin_quadric_tree_3<Kernel_>::Node::Operator::And));
}

template<class Kernel_>
bool Spin_quadric_tree_3<Kernel_>::evaluate(const Spin_3 &spin)
{
    assert(m_root);

    // FIXME: This is naive implementation
    return m_root->evaluate(spin);
}

template<class Kernel_>
size_t Spin_quadric_tree_3<Kernel_>::size() const
{
    if (m_root)
        return m_root->size();

    return 0;
}
} // namespace CS

#include "Spin_quadric_tree_3.ipp"

#endif // LIBCS_PREDICATE_SENTENCE_3_H
