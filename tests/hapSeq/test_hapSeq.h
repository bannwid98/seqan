// ==========================================================================
// %(TITLE)s
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: %(AUTHOR)s
// ==========================================================================

#ifndef SEQAN_TEST_HAPSEQ_H
#define SEQAN_TEST_HAPSEQ_H

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>

#include "../../apps/hapseq/variance.cpp"

SEQAN_DEFINE_TEST(test_isSNP_true)
{
    using namespace seqan2;

    VcfRecord record{};
    record.ref = "A";
    record.alt = "C";

    SEQAN_ASSERT(isSNP(record));

    record.alt = "C,G,T";

    SEQAN_ASSERT(isSNP(record));
}

SEQAN_DEFINE_TEST(test_isSNP_false)
{
    using namespace seqan2;

    VcfRecord record{};
    record.ref = "A";
    record.alt = "<DUP>";

    SEQAN_ASSERT(!isSNP(record));

    record.alt = "CGT";

    SEQAN_ASSERT(!isSNP(record));

    record.alt = "C,G,TA";

    SEQAN_ASSERT(!isSNP(record));

    record.ref = "ACGT";
    record.alt = "A";

    SEQAN_ASSERT(!isSNP(record));
}

SEQAN_DEFINE_TEST(test_isInDel_true)
{
    VcfRecord record{};
    record.alt = "ACGTN,AAAA,CCCC,TTTT,GGGG,NNNN";

    SEQAN_ASSERT(isInDel(record));
}

SEQAN_DEFINE_TEST(test_isInDel_false)
{
    VcfRecord record{};
    record.alt = "UUUU";

    SEQAN_ASSERT(!isInDel(record));

    record.alt = "<DUP>";

    SEQAN_ASSERT(!isInDel(record));

    record.alt = ".";

    SEQAN_ASSERT(!isInDel(record));

    record.alt = "";

    SEQAN_ASSERT(!isInDel(record));

    record.alt = "A[22:12530581[";

    SEQAN_ASSERT(!isInDel(record));
}

SEQAN_DEFINE_TEST(test_applySNP_valid)
{
    using namespace seqan2;

    Dna5String refSeq = "AAAAAAAAAAAAAAAAAAAAA";
    Dna5String hapSeq = "AAAAAAAAAAAAAAAAAAAAA";

    VcfRecord record{};
    record.beginPos = 10;
    record.ref = "A";
    record.alt = "C,G,T";

    unsigned altIdx = 1;
    hapSeq[10] = 'C';

    applySNP(refSeq, record, altIdx);
    SEQAN_ASSERT_EQ(refSeq, hapSeq);

    altIdx = 2;
    hapSeq[10] = 'G';

    applySNP(refSeq, record, altIdx);
    SEQAN_ASSERT_EQ(refSeq, hapSeq);

    altIdx = 3;
    hapSeq[10] = 'T';

    applySNP(refSeq, record, altIdx);
    SEQAN_ASSERT_EQ(refSeq, hapSeq);
}

SEQAN_DEFINE_TEST(test_applySNP_out_of_range) {
    Dna5String refSeq = "A";

    VcfRecord record{};
    record.beginPos = 10;
    record.ref = "A";
    record.alt = "U";

    bool caughtException{false};

    unsigned altIdx = 3;

    try {
        applySNP(refSeq, record, altIdx);
    } catch (std::out_of_range &e) {
        caughtException = true;
    }

    SEQAN_ASSERT(caughtException);

    caughtException = false;
    record.beginPos = 0;

    try {
        applySNP(refSeq, record, altIdx);
    } catch (std::out_of_range &e) {
        caughtException = true;
    }

    SEQAN_ASSERT(caughtException);
}

SEQAN_DEFINE_TEST(test_applySNP_invalid_argument) {
    using namespace seqan2;
    Dna5String refSeq = "A";

    VcfRecord record{};
    record.beginPos = 0;
    record.ref = "A";
    record.alt = "U";

    bool caughtException{false};

    unsigned altIdx = 1;

    try {
        applySNP(refSeq, record, altIdx);
    } catch (std::invalid_argument &e) {
        caughtException = true;
    }

    SEQAN_ASSERT(caughtException);
}

SEQAN_DEFINE_TEST(test_getAltIdx_valid)
{
    using namespace seqan2;

    VcfRecord record{};
    appendValue(record.genotypeInfos, "0");
    appendValue(record.genotypeInfos, "1");
    appendValue(record.genotypeInfos, "2");

    unsigned haploIdx = 0;
    unsigned altIdx = 0;

    SEQAN_ASSERT_EQ(getAltIdx(record, haploIdx), altIdx);

    haploIdx = 1;
    altIdx = 1;

    SEQAN_ASSERT_EQ(getAltIdx(record, haploIdx), altIdx);

    haploIdx = 2;
    altIdx = 2;

    SEQAN_ASSERT_EQ(getAltIdx(record, haploIdx), altIdx);
}

SEQAN_DEFINE_TEST(test_getAltIdx_out_of_range)
{
    VcfRecord record;
    appendValue(record.genotypeInfos, "1");
    appendValue(record.genotypeInfos, "2");

    bool exceptionThrown = false;

    try {
        getAltIdx(record, 3);
    } catch (std::out_of_range const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);
}

SEQAN_DEFINE_TEST(test_getInfo_valid) {
    VcfRecord record{};
    record.info = "SVTYPE=DUP;SVLEN=504;END=505;TARGETPOS=1:1328";

    std::unordered_map<std::string, std::string> infoMap{};

    getInfo(infoMap, record);

    SEQAN_ASSERT_EQ(infoMap["SVTYPE"], "DUP");
    SEQAN_ASSERT_EQ(infoMap["SVLEN"], "504");
    SEQAN_ASSERT_EQ(infoMap["END"], "505");
    SEQAN_ASSERT_EQ(infoMap["TARGETPOS"], "1:1328");

    record.info = "SVTYPE=DUP;SVTYPE=DEL;SVLEN=504;SVLEN=222;END=505;END=910;TARGETPOS=1:1328;TARGETPOS=8:123";
    infoMap.clear();

    std::streambuf* originalBuffer = std::cout.rdbuf();
    std::ostringstream oss;
    std::cout.rdbuf(oss.rdbuf());

    getInfo(infoMap, record);

    std::cout.rdbuf(originalBuffer);

    SEQAN_ASSERT(!oss.str().empty());

    SEQAN_ASSERT_EQ(infoMap["SVTYPE"], "DUP");
    SEQAN_ASSERT_EQ(infoMap["SVLEN"], "504");
    SEQAN_ASSERT_EQ(infoMap["END"], "505");
    SEQAN_ASSERT_EQ(infoMap["TARGETPOS"], "1:1328");
}

SEQAN_DEFINE_TEST(test_getInfo_parser_error) {
    VcfRecord record{};
    record.info = "SVTYPE=;SVLEN=;END=;TARGETPOS=;";

    std::unordered_map<std::string, std::string> infoMap{};

    bool exceptionThrown = false;

    try {
        getInfo(infoMap, record);
    } catch (ParseError const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);

    record.info = "SVTYPE:DUP;SVLEN:504;END:505;TARGETPOS:1:1328";
    exceptionThrown = false;

    try {
        getInfo(infoMap, record);
    } catch (ParseError const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);

    record.info = "SVTYPE=DUP,SVLEN=504,END=505,TARGETPOS=1:1328";
    exceptionThrown = false;

    try {
        getInfo(infoMap, record);
    } catch (ParseError const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);

    record.info = "SVTYPE=;SVLEN=;END=;TARGETPOS=";

    exceptionThrown = false;

    try {
        getInfo(infoMap, record);
    } catch (ParseError const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);
}

SEQAN_DEFINE_TEST(test_getInfo_empty) {
    VcfRecord record{};

    std::unordered_map<std::string, std::string> infoMap{};

    getInfo(infoMap, record);

    SEQAN_ASSERT(infoMap.empty());

    record.info = ".";

    getInfo(infoMap, record);

    SEQAN_ASSERT(infoMap.empty());

    record.info = "";

    getInfo(infoMap, record);

    SEQAN_ASSERT(infoMap.empty());
}

SEQAN_DEFINE_TEST(test_getEndPos_valid)
{
    std::unordered_map<std::string, std::string> infoMap = {{"END", "200"}, {"SVLEN", "100"}};
    size_t beginPos = 99;

    size_t endPos = 200;

    SEQAN_ASSERT_EQ(getEndPos(beginPos, infoMap), endPos);

    infoMap.clear();
    infoMap = {{"SVLEN", "100"}};

    SEQAN_ASSERT_EQ(getEndPos(beginPos, infoMap), endPos);

    infoMap.clear();
    infoMap =  {{"END", "200"}};

    SEQAN_ASSERT_EQ(getEndPos(beginPos, infoMap), endPos);
}

SEQAN_DEFINE_TEST(test_getEndPos_parser_error) {
    std::unordered_map<std::string, std::string> infoMap{{"OTHER", "100"}};
    size_t beginPos = 999999;

    bool exceptionThrown = false;

    try {
        getEndPos(beginPos, infoMap);
    } catch (ParseError const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);

    infoMap.clear();

    exceptionThrown = false;

    try {
        getEndPos(beginPos, infoMap);
    } catch (ParseError const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);
}

SEQAN_DEFINE_TEST(test_getEndPos_out_of_range) {
    std::unordered_map<std::string, std::string> infoMap{{"END", "0"}, {"SVLEN", "0"}};
    size_t beginPos = 99999;

    bool exceptionThrown = false;

    try {
        getEndPos(beginPos, infoMap);
    } catch (std::out_of_range const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);
    }

SEQAN_DEFINE_TEST(test_getEndPos_invalid_argument) {
    std::unordered_map<std::string, std::string> infoMap{{"END", "invalid"}};
    size_t beginPos = 100;

    bool exceptionThrown = false;

    try {
        getEndPos(beginPos, infoMap);
    } catch (std::invalid_argument const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);

    infoMap.clear();
    infoMap = {{"SVLEN", "invalid"}};

    exceptionThrown = false;

    try {
        getEndPos(beginPos, infoMap);
    } catch (std::invalid_argument const & e) {
        exceptionThrown = true;
    }

    SEQAN_ASSERT(exceptionThrown);
}

#endif
