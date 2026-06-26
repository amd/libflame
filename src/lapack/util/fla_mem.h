/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#ifndef FLA_MEM_H
#define FLA_MEM_H

#include "FLAME.h"

#define FLA_ALIGN(n, alignment) (((n) + (alignment)-1) & ~((alignment)-1))

/**
 * Allocates memory with the specified alignment.
 *
 * @param size The size of the memory block to allocate in bytes.
 * @param alignment The alignment in bytes (must be a power of two)
 * @return A pointer to the allocated memory block, or NULL if the allocation fails.
 */
void *fla_aligned_malloc(size_t size, size_t alignment);

/**
 * Frees memory allocated by fla_aligned_malloc.
 *
 * @param ptr A pointer to the memory block to free. This should be a pointer
 *            returned by a previous call to fla_aligned_malloc.
 */
void fla_aligned_free(void *ptr);

#endif /* FLA_MEM_H */