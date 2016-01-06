#ifndef __UPPSET_STOP_WATCH_HPP__
#define __UPPSET_STOP_WATCH_HPP__

#include <time.h>
#include <vector>
#include <sys/time.h>

namespace metag {

#define TIME_REPORT(rank, msg) { MPI_Barrier(MPI_COMM_WORLD); \
    double tm = clock.stop(); if (rank() == 0) { \
    cout << msg << ": " << tm << endl; } \
    clock.reset(); clock.start(); sync(); }


//--------------------------------------------------------------
//
//! \brief A generally useful class for measuring wallclock time
//
class StopWatch {
private:
    std::vector<struct timeval> m_lap;

public:
    StopWatch(size_t expectedLaps=1)  {
        reset( expectedLaps );
    }

    void reset(size_t expectedLaps=1) {
        m_lap.clear();
        m_lap.reserve( 2*expectedLaps );
    }

    bool start() {
        if (isRunning()) {
            return false;
        }
        struct timeval temp;
        if ( ::gettimeofday( &temp, 0 ) ) {
            return false;
        }
        m_lap.push_back( temp );
        return true;
    }

    bool lap() {
        struct timeval temp;
        if ( ::gettimeofday(&temp,0) ) {
            return false;
        }
        m_lap.push_back( temp );
        if (!isRunning()) {
            m_lap.push_back( temp );
        }
        return true;
    }
    double stop() {
        if (!isRunning()) {
            return -1.0;
        }
        struct timeval temp;
        if ( ::gettimeofday(&temp,0) ) {
            return -1.0;
        }
        m_lap.push_back( temp );
        return queryElapsedTime( );
    }

    double queryElapsedTime( int startLap=0, int lastLap=-1 ) const {
        int nlaps = queryLaps();
        if ( lastLap < 0 ) { // if negative, count backwards from the end.
            lastLap = nlaps + lastLap + 1;
        }
        if ( startLap < 0 ) { // if negative, count backwards from the end.
            startLap = nlaps + startLap + 1;
        }

        // now accumulate time intervals
        double time = 0.0;
        for ( int i=startLap; i<lastLap; ++i ) {
            double startsecs = m_lap[2*i  ].tv_sec + (1e-6 * m_lap[2*i  ].tv_usec);
            double endsecs =   m_lap[2*i+1].tv_sec + (1e-6 * m_lap[2*i+1].tv_usec);
            time += endsecs - startsecs;
        }
        return time;
    }

    bool isRunning() const {
        return (m_lap.size() % 2) == 1;
    }

    int queryLaps() const {
        return m_lap.size()/2;
    }

    double queryTotalTime() const {
        return queryElapsedTime();
    }

};

}
#endif
