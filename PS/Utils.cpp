#include "Utils.hpp"
#include <ostream>
#include <iostream>
#include <ctime>
#include <cstring>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

using std::cin;
using std::cout;
using std::string;

const int Bufsize = 1024;

void WaitForEnter()
{
    char c;
    while ( cin.rdbuf()->in_avail() ) cin.get(c);
    cout << "Press Enter to continue ...";
    c = '\0';
    while ( c != '\n' ) cin.get(c);
}


string GetFilename()
{
    char c;
    while ( cin.rdbuf()->in_avail() ) cin.get(c);
    cout << "Input file: ";
    string s;
    getline(cin,s);
    return s;
}


bool GetFilename(char* filename)
{
    char c;
    while ( cin.rdbuf()->in_avail() ) cin.get(c);
    cout << "Input file: ";
    cin.getline(filename,Bufsize);
    if ( strlen(filename)==0 ) {
        cout << "No file name specified \n";
        return false;
    }
    return true;
}


void Timer::Start()
{
    start = std::chrono::system_clock::now();
}

void Timer::Stop()
{
    end = std::chrono::system_clock::now();
}

void Timer::Display()
{
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    auto total_duration = end - start;
    auto seconds =
        std::chrono::duration_cast<std::chrono::seconds>(total_duration);
    auto milli =
        std::chrono::duration_cast<std::chrono::milliseconds>(total_duration
                                                              - seconds);
    Out << "Finished computation at " << std::ctime(&end_time);
    Out << "Elapsed time: " << seconds.count()
        << "s " << milli.count() << "ms\n";
}


int DisplayedCount::operator++()
{
    ++Count;
    if ( Count%Interval==0 ) {
        Out << Count;
        if (Count < 10*Interval) Out << "  ";
        else Out << " ";
        if (Count > 0 && Count%(10*Interval)==0) Out << '\n';
        Out.flush();
    }
    return Count;
}

int DisplayedCount::Final() const
{
    if (Count%Interval!=0) Out << Count << '\n';
    Out << '\n';
    return Count;
}


ProgressBar::ProgressBar(std::ostream& ostr, int total)
         : Out(ostr), Total(total), Interval(Total/100), Count(0)
{
    if ( Interval == 0 ) Interval = 1;
    Out << "0%       10        20        30        40        50        60        70        80        90        100%\n";
    Out << "|--------|---------|---------|---------|---------|---------|---------|---------|---------|---------|\n";
}

int ProgressBar::operator++()
{
    ++Count;
    if (Count%Interval==0) {
        Out << '*';
        Out.flush();
    }
    return Count;
}

int ProgressBar::Final() const
{
    if (Count%Interval!=0) Out << '*';
    Out << '\n';
    return Count;
}
