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

#define HASH_KEY(k) k // (k + (k >> 1))

template<class V, int P = int(1e7) + 19> struct Hash {
    ListNode<V> *ptr[P];
    int count = 0;

    void clear() {
        // buf.reset();
        count = 0;
        for (int i = 0; i < P; ++i) {
            ptr[i] = nullptr;
        }
    }

    int validEntries() {
        int x = 0;
        for (int i = 0; i < P; ++i) {
            x += !!ptr[i];
        }
        return x;
    }

    ListNode<V> *find(unsigned long key) {
        return ptr[HASH_KEY(key) % P];
    }

    void insert(unsigned long key, ListNode<V> *node) {
        // if (key == 4385068548080ll) {
        //     printf("yes!\n");
        // }

        int loc = HASH_KEY(key) % P;
        node->next = ptr[loc];
        ptr[loc] = node;
    }
};
