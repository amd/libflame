/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "fla_mem.h"

/**
 * Allocates memory with the specified alignment.
 *
 * @param size The size of the memory block to allocate in bytes.
 * @param alignment The alignment in bytes (must be a power of two and at least the size of a
 * pointer)
 * @return A pointer to the allocated memory block, or NULL if the allocation fails.
 */
void *fla_aligned_malloc(size_t size, size_t alignment)
{
    if(alignment < sizeof(void *) || (alignment & (alignment - 1)))
    {
        // Alignment must be a power of two and at least the size of a pointer
        return NULL;
    }
    size_t align_mask = alignment - 1;

    size_t padding = align_mask + sizeof(void *);

    void *raw_ptr = malloc(size + padding);
    if(raw_ptr == NULL)
    {
        return NULL;
    }

    uintptr_t raw_addr = (uintptr_t)raw_ptr;

    uintptr_t aligned_addr = (raw_addr + padding) & ~align_mask;
    void **aligned_ptr = (void **)aligned_addr;
    aligned_ptr[-1] = raw_ptr;
    return (void *)aligned_ptr;
}

/**
 * Frees memory allocated by fla_aligned_malloc.
 *
 * @param ptr A pointer to the memory block to free. This should be a pointer
 *            returned by a previous call to fla_aligned_malloc.
 */
void fla_aligned_free(void *ptr)
{
    if(ptr == NULL)
    {
        return;
    }
    void *raw_ptr = ((void **)ptr)[-1];
    free(raw_ptr);
}