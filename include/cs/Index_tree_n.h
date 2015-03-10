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
#ifndef LIBCS_INDEX_TREE_N_H
#define LIBCS_INDEX_TREE_N_H

#include <log4cxx/logger.h>
#include <cstdint>
#include <cassert>

namespace CS
{
struct Empty_node_data {};

template<class Cell_, class Node_data_>
class Index_tree_node_n
{
    typedef typename Cell_::Handle Cell_handle;

public:
    class Pool
    {
        typedef Index_tree_node_n<Cell_, Node_data_> Index_tree_node;

    public:
        // Note: this may be too few in some extreme cases!
        typedef boost::uint32_t     Offset;

        // null offset which does not point anything
        static const Offset NULL_OFFSET = static_cast<Offset>(-1);

        Pool(size_t total)
            : m_logger(log4cxx::Logger::getLogger("CS.Index_tree_node_n.Pool"))
        {
            LOG4CXX_DEBUG(m_logger, "Creating a pool for " << total << " structures of size " << sizeof(Index_tree_node) << " bytes");
            LOG4CXX_DEBUG(m_logger, "Maximum pool size is " << std::fixed << std::setprecision(2) << ((total *  sizeof(Index_tree_node)) / (1024 * 1024)) << " MB");

            // the pool is very mean at memory allocation - it does allocate only
            // when needed - always increase the side by a factor of two
            m_allocated = total >> 3; // take 1/8 of all requested memory

            if (m_allocated < 1)
                m_allocated = 1;

            size_t data_size = m_allocated * sizeof(Index_tree_node);

            m_data = static_cast<Index_tree_node *>(malloc(data_size));
            assert(m_data);

            m_total = total;
            m_used = 0;

            LOG4CXX_DEBUG(m_logger, "New allocated pool size is " << std::fixed << std::setprecision(2) << (data_size / (1024 * 1024)) << " MB");
        }

        ~Pool()
        {
            // kill all memory
            free(m_data);
        }

        // allocate new pool index
        Offset alloc()
        {
            // never allow more allocations than declared
            assert(m_used < m_total);

            // adjust allocated size
            if (m_used == m_allocated)
            {
                // double the size
                m_allocated <<= 1;

                // clamp to maximum
                if (m_allocated > m_total)
                    m_allocated = m_total;

                // reallocate
                size_t data_size = m_allocated * sizeof(Index_tree_node);

                LOG4CXX_DEBUG(m_logger, "Reallocating pool...");

                Index_tree_node *data = static_cast<Index_tree_node *>(realloc(m_data, data_size));
                assert(data);

                m_data = data;

                LOG4CXX_DEBUG(m_logger, "New allocated pool size is " << std::fixed << std::setprecision(2) << (data_size / (1024 * 1024)) << " MB");
            }

            return m_used++;
        }

        Index_tree_node *ref(Offset offset)
        {
            assert(offset != NULL_OFFSET);
            return m_data + offset;
        }

        const Index_tree_node *ref(Offset offset) const
        {
            assert(offset != NULL_OFFSET);
            return m_data + offset;
        }

        double usage() const
        {
            return 100.0 * m_used / m_total;
        }

        size_t memsize() const
        {
            return m_total * sizeof(Index_tree_node);
        }

        size_t allocated_memsize() const
        {
            return m_allocated * sizeof(Index_tree_node);
        }

    private:
        Index_tree_node *   m_data;
        size_t              m_total;
        size_t              m_used;
        size_t              m_allocated;

        log4cxx::LoggerPtr  m_logger;
    };

    // an offset to the node in an integrated pool allocator
    typedef typename Pool::Offset Offset;

    // null offset which does not point anything
    static const Offset NULL_OFFSET = Pool::NULL_OFFSET;

    static Offset create_root(Pool *pool)
    {
        // WARN: this must be a separate call because of the posible reallocation
        Offset root = pool->alloc();

        pool->ref(root)->construct(NULL_OFFSET);
        return root;
    }

    template<typename InputIterator>
    static Offset store(InputIterator begin, InputIterator end, Cell_handle cell_handle, Pool *pool, Offset root)
    {
        Offset current = root;
        int bit;

        while (begin != end)
        {
            bit = !!*begin;

            if (pool->ref(current)->m_child[bit] == NULL_OFFSET)
            {
                // WARN: this must be a separate call because of the posible reallocation
                Offset new_offset = pool->alloc();

                pool->ref(current)->m_child[bit] = new_offset;
                pool->ref(pool->ref(current)->m_child[bit])->construct(current);
            }

            current = pool->ref(current)->m_child[bit];
            ++begin;
        }

        // store cell coodinate
        pool->ref(current)->m_cell_handle = cell_handle;
        return current;
    }

    template<typename InputIterator>
    static std::pair<Offset, bool> check_and_store(InputIterator begin, InputIterator end, Pool *pool, Offset root)
    {
        InputIterator next = begin;
        Offset node = root;
        bool flag = false; // new node created flag
        int bit;

        while (begin != end)
        {
            ++next; // next node is now accessible

            bit = !!*begin;

            if (pool->ref(node)->m_child[bit] == NULL_OFFSET)
            {
                if (next == end)
                    flag = true;

                // WARN: this must be a separate call because of the posible reallocation
                Offset new_offset = pool->alloc();

                pool->ref(node)->m_child[bit] = new_offset;
                pool->ref(pool->ref(node)->m_child[bit])->construct(node);
            }

            node = pool->ref(node)->m_child[bit];
            ++begin;
        }

        return std::make_pair(node, flag);
    }

    template<typename InputIterator>
    static boost::optional<Cell_handle> contains(InputIterator begin, InputIterator end, Pool *pool, Offset root)
    {
        if (begin == end)
            return pool->ref(root)->m_cell_handle;

        int bit = !!*begin;

        if (pool->ref(root)->m_child[bit] == NULL_OFFSET)
            return boost::optional<Cell_handle>();

        // TODO: unwind stack
        return contains(++begin, end, pool, pool->ref(root)->m_child[bit]);
    }

    static bool is_left(Pool *pool, Offset node)
    {
        return !pool->ref(node)->is_root() &&
                pool->ref(pool->ref(node)->parent())->left() == node;
    }

    static bool is_right(Pool *pool, Offset node)
    {
        return !pool->ref(node)->is_root() &&
                pool->ref(pool->ref(node)->parent())->right() == node;
    }

    bool is_root() const
    {
        return m_parent == NULL_OFFSET;
    }

    bool has_left() const
    {
        return left() != NULL_OFFSET;
    }

    bool has_right() const
    {
        return right() != NULL_OFFSET;
    }

    Offset left() const
    {
        return m_child[0];
    }

    Offset right() const
    {
        return m_child[1];
    }

    Offset parent() const
    {
        return m_parent;
    }

    bool is_leaf() const
    {
        return !has_left() && !has_right();
    }

    Cell_handle cell_handle() const
    {
        return m_cell_handle;
    }

    void set_cell_handle(Cell_handle cell_handle)
    {
        m_cell_handle = cell_handle;
    }

    Node_data_ data() const
    {
        return m_data;
    }

    void set_data(Node_data_ data)
    {
        m_data = data;
    }

private:
    Index_tree_node_n() {}

    void construct(Offset parent)
    {
        m_child[0] = NULL_OFFSET;
        m_child[1] = NULL_OFFSET;
        m_parent = parent;
        m_data = Node_data_();
    }

    Offset      m_child[2];
    Offset      m_parent;
    Cell_handle m_cell_handle;
    Node_data_   m_data;
};

template<class Cell_>
struct Index_tree_node_visitor_n
{
    typedef typename Cell_::Handle Cell_handle;
    typedef Index_tree_node_n<Cell_, Index_tree_node_visitor_n *> Index_tree_node;

    Index_tree_node_visitor_n(
            typename Index_tree_node::Offset node_,
            Cell_handle handle_)
        : node(node_),
          handle(handle_)
    {
    }

    // current level node
    typename Index_tree_node::Offset node;

    // final cell handle
    Cell_handle handle;
};

template<class Value_>
class Multi_list_link_n
{
public:
    const Value_ &       value() const
    {
        return m_value;
    }

    Value_ &             value()
    {
        return m_value;
    }

    void                set_value(const Value_ &value)
    {
        m_value = value;
    }

    Multi_list_link_n * prev() const
    {
        return m_prev;
    }

    Multi_list_link_n * next() const
    {
        return m_next;
    }

private:
    Value_               m_value;
    Multi_list_link_n * m_prev;
    Multi_list_link_n * m_next;

    template<class OtherValue>
    friend class Multi_list_n;

    template<class OtherValue>
    friend class Multi_list_rope_n;
};

template<class Value_>
class Multi_list_rope_n
{
    typedef Multi_list_link_n<Value_> Multi_list_link;

public:
    void                push_link(Multi_list_link *link)
    {
        link->m_next = m_links_head;
        if (m_links_head) m_links_head->m_prev = link;
        link->m_prev = 0;
        m_links_head = link;
    }

    void                pop_link(Multi_list_link *link)
    {
        if (link->m_next) link->m_next->m_prev = link->m_prev;
        if (link->m_prev) link->m_prev->m_next = link->m_next;
        if (m_links_head == link) m_links_head = link->m_next;
    }

    Multi_list_link *   links_head() const
    {
        return m_links_head;
    }

    Multi_list_rope_n * prev() const
    {
        return m_prev;
    }

    Multi_list_rope_n * next() const
    {
        return m_next;
    }

private:
    Multi_list_link *   m_links_head;

    Multi_list_rope_n * m_prev;
    Multi_list_rope_n * m_next;

    template<class OtherValue>
    friend class Multi_list_n;
};

template<class Value_>
class Multi_list_n
{
    typedef Multi_list_link_n<Value_> Multi_list_link;
    typedef Multi_list_rope_n<Value_> Multi_list_rope;

public:
    Multi_list_n(size_t size)
    {
        m_size = size;

        m_rope_pool = static_cast<Multi_list_rope *>(malloc(m_size * sizeof(Multi_list_rope)));
        m_link_pool = static_cast<Multi_list_link *>(malloc(m_size * sizeof(Multi_list_link)));

        m_free_ropes = 0;
        m_free_links = 0;

        for (size_t i = 0; i < m_size; ++i)
        {
            m_rope_pool[i].m_next = m_free_ropes;
            m_free_ropes = &m_rope_pool[i];

            m_link_pool[i].m_next = m_free_links;
            m_free_links = &m_link_pool[i];
        }

        m_ropes_head = 0;
    }

    ~Multi_list_n()
    {
        free(m_rope_pool);
        free(m_link_pool);
    }

    size_t memsize() const
    {
        return m_size * (sizeof(Multi_list_rope) + sizeof(Multi_list_link));
    }

    Multi_list_link *alloc_link()
    {
        assert(m_free_links);
        Multi_list_link *result = m_free_links;
        m_free_links = m_free_links->m_next;
        result->m_prev = result->m_next = 0;
        return result;
    }

    void release_link(Multi_list_link *link)
    {
        link->m_next = m_free_links;
        m_free_links = link;
    }

    Multi_list_rope *alloc_rope()
    {
        assert(m_free_ropes);
        Multi_list_rope *result = m_free_ropes;
        m_free_ropes = m_free_ropes->m_next;
        result->m_prev = result->m_next = 0;
        result->m_links_head = 0;
        return result;
    }

    void release_rope(Multi_list_rope *rope)
    {
        rope->m_next = m_free_ropes;
        m_free_ropes = rope;
    }

    void push_rope(Multi_list_rope *rope)
    {
        rope->m_next = m_ropes_head;
        if (m_ropes_head) m_ropes_head->m_prev = rope;
        rope->m_prev = 0;
        m_ropes_head = rope;
    }

    void pop_rope(Multi_list_rope *rope)
    {
        if (rope->m_next) rope->m_next->m_prev = rope->m_prev;
        if (rope->m_prev) rope->m_prev->m_next = rope->m_next;
        if (m_ropes_head == rope) m_ropes_head = rope->m_next;
    }

    Multi_list_rope *ropes_head() const
    {
        return m_ropes_head;
    }

private:
    // total size of ropes and links
    size_t      m_size;

    // pools
    Multi_list_rope *       m_rope_pool;
    Multi_list_link *       m_link_pool;

    // free
    Multi_list_rope *       m_free_ropes;
    Multi_list_link *       m_free_links;

    // ropes
    Multi_list_rope *       m_ropes_head;
};
} // namespace CS

#include "Index_tree_n.ipp"

#endif // LIBCS_INDEX_TREE_N_H
