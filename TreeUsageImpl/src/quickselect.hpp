#ifndef QUICKSELECT_HPP
#define QUICKSELECT_HPP

#include <iostream>
#include <string>
#include <vector>


void swap(double *const a, double *const b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

unsigned long long partition(std::vector<double> &arr, const unsigned long long low, const unsigned long long high)
{
    const double pivot = arr[high];
    unsigned long long i = low - 1;

    for (unsigned long long j = low; j < high; j++)
    {
        if (arr[j] <= pivot)
        {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }

    swap(&arr[i + 1], &arr[high]);
    return i + 1;
}

double quickSelect(std::vector<double> &arr, const unsigned long long low, const unsigned long long high, const int k)
{
    if (low <= high)
    {
        const int pivotIndex = partition(arr, low, high);

        if (pivotIndex == k)
        {
            return arr[pivotIndex];
        }
        else if (pivotIndex < k)
        {
            return quickSelect(arr, pivotIndex + 1, high, k);
        }
        else
        {
            return quickSelect(arr, low, pivotIndex - 1, k);
        }
    }

    return -1.0;
}


#endif