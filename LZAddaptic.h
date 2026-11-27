/////////////////////////////////////////////////////////////////
//
//	LZAddaptic.h
//
//	Part of Addaptic
//
//	Copyright (c) 2000-2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

#pragma once


/////////////////////////////////////////////////////////////////
//  All versions
/////////////////////////////////////////////////////////////////

#include <cstdlib>		// calloc, free
#include <cstring>		// the string control macros below


/////////////////////////////////////////////////////////////////
//  String control
/////////////////////////////////////////////////////////////////

#define STRINGLENGTH ::strlen
#define STRINGCOPY ::strcpy
#define STRINGCMP ::strcmp


/////////////////////////////////////////////////////////////////
//  Addaptic version constants
//
//  Numbered by semantic versioning: the major number moves for a
//  change that can break a caller, the minor for something added,
//  and the release for a fix. A fix may change an answer -- one that
//  was not correctly rounded becoming one that is -- and is still a
//  fix. CMakeLists.txt gives the same number as the project's
//  version, and stops if the two disagree.
/////////////////////////////////////////////////////////////////

constexpr long		kAddapticVersionMajor = 1;
constexpr long		kAddapticVersionMinor = 0;
constexpr long		kAddapticVersionRelease = 0;

constexpr char		kAddapticVersionString[] = "1.0.0";


/////////////////////////////////////////////////////////////////
//  Error codes
/////////////////////////////////////////////////////////////////

const long	kNoError = 0;
