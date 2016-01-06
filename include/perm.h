#ifndef __PERM_H__
#define __PERM_H__

#if WITH_PJMALLOC == 1
#define JEMALLOC_MANGLE // use if configured --with-jemalloc-prefix=pj
// #define PERM_OVERRIDE_NEW // use to override new and delete operators




#include "jemalloc/pallocator.h"


class Arena {
protected:
	void *m_heap;

public:

	template <typename T>
	T* create()
		{ return PERM_NEW(T)(); }

	template <typename T, typename A0>
	T* create(const A0& _a0)
		{ return PERM_NEW(T)(_a0); }

	template <typename T, typename A0, typename A1>
	T* create(const A0& _a0, const A1& _a1)
		{ return PERM_NEW(T)(_a0, _a1); }

	template <typename T, typename A0, typename A1, typename A2>
	T* create(const A0& _a0, const A1& _a1, const A2& _a2)
		{ return PERM_NEW(T)(_a0, _a1, _a2); }

	template <typename T>
	void destroy(T* _addr)
		{ PERM_DELETE(_addr,T); }

};

class Persistent : public Arena { // derived from Arena
private:
	/* globals */

public:

	/* Register a block as persistent memory */
	int perm(void *ptr, size_t size)
	 { return ::perm(ptr, size); }

	/* Open and map file into core memory */
	int mopen(const char *fname, const char *mode, size_t size)
		{ m_heap = NULL; return ::mopen(fname, mode, size); }

	/* Close memory-mapped file */
	int mclose(void)
		{ return ::mclose(); }

	/* Flushes in-core data to memory-mapped file */
	int mflush(void)
		{ return ::mflush(); }

	/* Open backup file */
	int bopen(const char *fname, const char *mode);

	/* Close backup file */
	int bclose(void);

	/* Backup globals and heap to a separate file */
	int backup(void);

	/* Restore globals and heap from a separate file */
	int restore(void);

};

#else // WITH_PJMALLOC == 1

class Arena {};
class Persistent : public Arena {};

#endif // WITH_PJMALLOC == 1

#endif
