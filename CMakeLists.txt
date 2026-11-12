cmake_minimum_required(VERSION 3.24)
project(Addaptic)

set(CMAKE_CXX_STANDARD 17)

# The suite is an accuracy gate and will not notice a stray read or
# write into the digits array - those land in the slack the array is
# calloc'd with, and change nothing observable. Configure with
# -DADDAPTIC_ASAN=ON to put a memory-safety gate behind it.
option(ADDAPTIC_ASAN "Build with AddressSanitizer" OFF)

# The threads gate sees a race that changes an answer, crashes or hangs.
# One that leaves every answer right only ThreadSanitizer sees, so
# configure with -DADDAPTIC_TSAN=ON after anything that touches state a
# calculation keeps. It cannot be combined with AddressSanitizer.
option(ADDAPTIC_TSAN "Build with ThreadSanitizer" OFF)

if(ADDAPTIC_ASAN AND ADDAPTIC_TSAN)
    message(FATAL_ERROR "ADDAPTIC_ASAN and ADDAPTIC_TSAN cannot be combined")
endif()

set(ADDAPTIC_ENGINE
        LZMath.cpp
        LZNumber.cpp)

# MSVC spells warnings its own way, and its -Wall is every warning the
# system headers carry as well. /W4 is its nearest to -Wall -Wextra, and
# /permissive- holds it to the standard, so code it accepts any MSVC
# configuration will -- Elzed's included.
if(MSVC)
    set(ADDAPTIC_WARNINGS /W4 /permissive-)
else()
    set(ADDAPTIC_WARNINGS -Wall -Wextra)
endif()

# MSVC's AddressSanitizer is a compiler switch alone: its objects name the
# runtime for the linker themselves, and link.exe only warns at
# /fsanitize. It cannot share a build with /RTC1, the run-time checks
# CMake gives Debug, so they go, and /Z7 keeps the debug information in
# the objects, so targets compiling the same source never contend for a
# PDB. To run, clang_rt.asan_dynamic-x86_64.dll has to be on the PATH --
# it sits beside cl.exe. MSVC has no ThreadSanitizer.
if(ADDAPTIC_ASAN)
    if(MSVC)
        set(ADDAPTIC_SANITIZE /fsanitize=address /Z7)
        set(ADDAPTIC_SANITIZE_LINK /INCREMENTAL:NO /DEBUG)
        string(REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
        string(REPLACE "/Zi" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    else()
        set(ADDAPTIC_SANITIZE -fsanitize=address -fno-omit-frame-pointer -g)
        set(ADDAPTIC_SANITIZE_LINK ${ADDAPTIC_SANITIZE})
    endif()
endif()

if(ADDAPTIC_TSAN)
    if(MSVC)
        message(FATAL_ERROR "MSVC has no ThreadSanitizer")
    endif()
    set(ADDAPTIC_SANITIZE -fsanitize=thread -g)
    set(ADDAPTIC_SANITIZE_LINK ${ADDAPTIC_SANITIZE})
endif()

add_executable(Addaptic
        Addaptic.cpp
        ${ADDAPTIC_ENGINE})

target_compile_options(Addaptic PRIVATE ${ADDAPTIC_WARNINGS} ${ADDAPTIC_SANITIZE})
target_link_options(Addaptic PRIVATE ${ADDAPTIC_SANITIZE_LINK})

add_executable(AddapticTests
        AddapticTests.cpp
        ${ADDAPTIC_ENGINE})

target_compile_options(AddapticTests PRIVATE ${ADDAPTIC_WARNINGS} ${ADDAPTIC_SANITIZE})
target_link_options(AddapticTests PRIVATE ${ADDAPTIC_SANITIZE_LINK})

# The magnitude sweep. AddapticTests measures how good an answer is;
# this asks whether every function still answers as well as it did at
# every size, from 1e-40 to 1e+40. Seven silent wrong answers were
# found in one week that the accuracy suite could not see, because it
# had nothing above 1e30 or below 1e-7.
add_executable(AddapticSweep
        AddapticSweep.cpp
        ${ADDAPTIC_ENGINE})

target_compile_options(AddapticSweep PRIVATE ${ADDAPTIC_WARNINGS} ${ADDAPTIC_SANITIZE})
target_link_options(AddapticSweep PRIVATE ${ADDAPTIC_SANITIZE_LINK})

# The threads gate. Every census function computed by eight threads at
# once, each answer compared with the one a single thread gives, with
# each thread's own operands and with operands two threads share; and a
# new thread's working precision checked.
find_package(Threads REQUIRED)

add_executable(AddapticThreads
        AddapticThreads.cpp
        ${ADDAPTIC_ENGINE})

target_compile_options(AddapticThreads PRIVATE ${ADDAPTIC_WARNINGS} ${ADDAPTIC_SANITIZE})
target_link_options(AddapticThreads PRIVATE ${ADDAPTIC_SANITIZE_LINK})
target_link_libraries(AddapticThreads PRIVATE Threads::Threads)

# The audit's binary. auditrounding.py hands it the sweep's cases and
# scores what comes back; it is built with the rest so that every
# compiler that builds the engine builds this too -- the Windows run
# compares its raw answers with the Mac's.
add_executable(auditrounding
        auditrounding.cpp
        ${ADDAPTIC_ENGINE})

target_compile_options(auditrounding PRIVATE ${ADDAPTIC_WARNINGS} ${ADDAPTIC_SANITIZE})
target_link_options(auditrounding PRIVATE ${ADDAPTIC_SANITIZE_LINK})

enable_testing()
add_test(NAME accuracy COMMAND AddapticTests)
add_test(NAME magnitudes COMMAND AddapticSweep)
add_test(NAME threads COMMAND AddapticThreads)

# ASan reports through the exit status, so the accuracy test doubles as
# the memory-safety test when the sanitizer is on...
if(ADDAPTIC_ASAN)
    set_tests_properties(accuracy magnitudes threads PROPERTIES
            ENVIRONMENT "ASAN_OPTIONS=abort_on_error=1:detect_leaks=0")
endif()

# A fault at the extremes has twice been a hang rather than a wrong
# answer, and a hung test that never returns is worse than a failing
# one. Neither of these should take more than a couple of seconds.
set_tests_properties(accuracy magnitudes PROPERTIES TIMEOUT 120)

# The threads gate takes a few seconds, and ThreadSanitizer slows it
# many times over...
set_tests_properties(threads PROPERTIES TIMEOUT 900)
