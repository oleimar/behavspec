#ifndef UTILS_HPP
#define UTILS_HPP

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

/***************************************************************************
Utils.hpp

This unit provides a number of often used functions and classes, mainly for
console applications, such as prompting for a file name, displaying a counter
or a progress bar, and displaying the time elapsed.

***************************************************************************/

#include <chrono>
#include <string>
#include <iosfwd>


// WaitForEnter: Pause until user presses enter
void WaitForEnter();


// GetFilename: Prompt user for file name and return it as a string
std::string GetFilename();


// GetFilename: Prompt user for file name; the buffer char* filename should
// have room for Bufsize characters (including terminating \0)
bool GetFilename(char* filename);


// Timer: A "stop watch" class; the stream on which to display
// the time (usually std::cout) is passed to the constructor;
class Timer {
    std::ostream& Out;
    std::chrono::time_point<std::chrono::system_clock> start, end;
public:
    explicit Timer(std::ostream& ostr) : Out(ostr) {}
    void Start();   // start timing
    void Stop();    // stop timing
    void Display(); // display elapsed time between start and stop
};


// DisplayedCount: Use this class for a counter that displays every tenth
// on a stream (usually std::cout)
class DisplayedCount {
    std::ostream& Out;
    int Count;
    int Interval;
public:
    explicit DisplayedCount(std::ostream& ostr, int count = 0, int interval = 10)
            : Out(ostr), Count(count), Interval(interval) { }
    // Increment and display if count is a multiple of interval
    int operator++();
    operator int() const { return Count; } // Conversion to int operator
    // Display final count (unless it is a multiple of interval)
    int Final() const;
};


// ProgressBar: Use this class for a progress bar that displays percentages
// from 0 to 100 on a stream (usually std::cout)
class ProgressBar {
    std::ostream& Out;
    int Total;
    int Interval;
    int Count;
public:
    explicit ProgressBar(std::ostream& ostr, int total);
    int operator++(); // Increment count and display if it is a multiple of total/100
    operator int() const { return Count; } // Conversion to int operator
    int Final() const; // Display final
};


#endif // UTILS_HPP
