#include <cassert>
using namespace std;

template<class V> struct ListNode {
    unsigned long key;
    V item;
    ListNode *next;

    ListNode() {}
};

template<class T, int BUFFER_SIZE = int(2e6)> struct Buffer {
    T buf[BUFFER_SIZE], *bufptr = buf;

    Buffer() {}

    T *get() {
        assert(bufptr - buf < BUFFER_SIZE);
        return bufptr++;
    }

    T *begin() {
        return buf;
    }

    T *end() {
        return bufptr;
    }

    void reset() {
        bufptr = buf;
    }

    void add(T item) {
        *bufptr = item;
        ++bufptr;
    }

    int size() {
        return bufptr - buf;
    }
};

template<class V, int P = int(1e7) + 19> struct Hash {
    ListNode<V> *ptr[P];

    void clear() {
        // buf.reset();
        for (int i = 0; i < P; ++i) {
            ptr[i] = nullptr;
        }
    }

    ListNode<V> *find(unsigned long key) {
        return ptr[key % P];
    }

    void insert(unsigned long key, ListNode<V> *node) {
        // if (key == 4385068548080ll) {
        //     printf("yes!\n");
        // }

        int loc = key % P;
        node->next = ptr[loc];
        ptr[loc] = node;
    }
};
