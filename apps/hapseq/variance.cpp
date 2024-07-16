//
// Created by dominik on 26.06.24.
//

#include <seqan/vcf_io.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <variant>
#include <filesystem>
#include <iostream>
#include <algorithm>

using namespace seqan2;

unsigned getAltIdx (const VcfRecord & record, const unsigned haploIdx) {

    if (haploIdx >= length(record.genotypeInfos)) {
        throw std::out_of_range("HaploIdx is out of range of genotypeInfos.");
    }

    return std::stoi(toCString(record.genotypeInfos[haploIdx]));
}

IsInAlphabet<Dna5> isInDNA5{};

struct ConsensusOptions {
    std::filesystem::path variancePath{};
    std::filesystem::path referencePath{};
    String<CharString> haplotypeNames{};
    std::filesystem::path outdirPath{};
};

void setupParser(ArgumentParser & parser) {
    setShortDescription(parser, "Generate hapseq sequences from a reference and a variance file.");
    setVersion(parser, "0.0.1");
    setDate(parser, "June 2024");

    addUsageLine(parser,
                 R"([\fIOPTIONS\fP] "\fITEXT\fP")");
    addDescription(parser,
                   "This program creates a hapseq sequences from a variance file in vcf format and reference"
                   "sequence in fasta format. It creates a sequence for each listed haplotype (-H) and saves them in "
                   "provided output directory as fasta files. Beside small SNPs and small indels also larger structural "
                   "variants, such as inversions, duplications and translocations are support. Since structural "
                   "variants can be complex and not well standardized in vcf format, please check the example files for "
                   "the supported format.");

    addSection(parser, "Mandatory Input");

    addOption(parser, ArgParseOption("v",
                                     "variance",
                                     "Path to input variance file in vcf format",
                                     ArgParseArgument::INPUT_FILE,
                                     "IN"));
    addOption(parser, ArgParseOption("r",
                                     "reference",
                                     "Path to input reference file in fasta format",
                                     ArgParseArgument::INPUT_FILE,
                                     "IN"));
    addOption(parser, ArgParseOption("H",
                                     "haplotype",
                                     "Name of haplotype to create hapseq sequence for. "
                                     "Can be given multiple times.",
                                     ArgParseArgument::STRING,
                                     "STR",
                                     true));

    addSection(parser, "Mandatory Output");

    addOption(parser, ArgParseOption("o",
                                     "outdir",
                                     "Path to directory, where the hapseq files get saved "
                                     "in fasta format. The files will have the name of the haplotypes.",
                                     ArgParseArgument::OUTPUT_DIRECTORY,
                                     "OUT"));

    addTextSection(parser, "Examples");

    addListItem(parser,
                R"(\fBconsensus\fP \fB-v\fP \fIpath/to/variance_file.vcf\fP \"
                "\fB-r\fP \fIpath/to/reference_file.fasta\fP \"
                "\fB-H\fP \fIhaplotype1\fP \fB-H\fP \fIhaplotype2\fP \"
                "\fB-o\fP \fIpath/to/output_directory\fP)",
                "Creates the sequence for haplotype1 and haplotype2 from the variance file and reference file "
                "and saves them in the output directory.");

    setRequired(parser, "outdir");
    setRequired(parser, "variance");
    setRequired(parser, "reference");
    setRequired(parser, "haplotype");

    setValidValues(parser, "variance", "vcf");
    setValidValues(parser, "reference", "fasta fa");
}

ArgumentParser::ParseResult getOptionsFromParser(const ArgumentParser & parser, ConsensusOptions & options) {

    getOptionValue(options.outdirPath, parser, "outdir");

    CharString haplotypeName{};

    for (size_t i = 0; i < getOptionValueCount(parser, "haplotype"); ++i) {
        getOptionValue(haplotypeName, parser, "haplotype", i);
        appendValue(options.haplotypeNames, haplotypeName);
    }

    getOptionValue(options.variancePath, parser, "variance");
    getOptionValue(options.referencePath, parser, "reference");

    return ArgumentParser::PARSE_OK;
}

ArgumentParser::ParseResult parseCommandLine(ConsensusOptions & options, int argc, char ** argv) {

    ArgumentParser parser("hapseq");

    setupParser(parser);

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ::ArgumentParser::PARSE_OK)
        return res;

    return getOptionsFromParser(parser, options);
}

bool loadOrBuildFAI(FaiIndex & faiIndex, const std::filesystem::path & path) {

    bool res = true;

    if (!open(faiIndex, path.c_str())) {

        if (!build(faiIndex, path.c_str())) {
            res = false;
        } else {
            CharString faiPath = path.c_str();
            append(faiPath, ".fai");
            try {
                save(faiIndex, toCString(faiPath));
            }
            catch (IOError const &ioErr) {
                res = false;
            }
        }
    }
    return res;
}

bool openAndLoadVCFHeader(VcfHeader & header, VcfFileIn & vcfIn, const std::filesystem::path & path) {

    bool res = true;

    if (!open(vcfIn, path.c_str())) {
        std::cerr << "Could not open variance file!" << std::endl;
        res = false;
    } else {
        try {
            readHeader(header, vcfIn);
        }
        catch (ParseError const &e) {
            std::cerr << "ERROR: input header is badly formatted. " << e.what() << "\n";
            res = false;
        }
        catch (IOError const &e) {
            std::cerr << "ERROR: could not read header. " << e.what() << "\n";
            res = false;
        }
    }

    return res;
}

bool isSNP(const VcfRecord & record) {

    bool res = true;

    bool prevWasComma = true;
    for (char c : record.alt) {
        if (c == ',') {
            prevWasComma = true;
        } else {
            if(!prevWasComma) {
                res = false;
                break;
            } else {
                prevWasComma = false;
            }
        }
    }

    return res && length(record.ref) == 1;
}

bool isInDel(const VcfRecord & record) {

    bool res = !empty(record.alt);

    for (char c : record.alt) {
        if (isInDNA5(c) || c == ',') {
            continue;
        } else {
            res = false;
            break;
        }
    }

    return res;
}

void getInfo(std::unordered_map<std::string, std::string> & infoMap, const VcfRecord & record) {
    std::istringstream iss{toCString(record.info)};
    std::string key, value;

    if (record.info == ".") {
        return;
    }

    while (std::getline(iss, key, '=')) {
        std::getline(iss, value, ';');

        if (key.empty() || value.empty() || std::find(std::begin(value), std::end(value), '=') != std::end(value)) {
            throw ParseError(
                    "ERROR: INFO dose not non empty contain key value pairs in form 'key=value' seperated by ';'.");
        }

        if (infoMap.find(key) != infoMap.end()) {
            std::cout << "WARNING: INFO field contains duplicate keys. Only the first one will be used." << std::endl;
        }

        infoMap.insert({key, value});
    }
}

void applySNP(Dna5String & seq, const VcfRecord & record, const unsigned altIdx) {

    if (record.beginPos >= length(seq)) {
        SEQAN_THROW(std::out_of_range("Position of variance is out of range of sequence"));
    }

    unsigned altPos{(altIdx-1)*2};

    if (altPos >= length(record.alt)) {
        SEQAN_THROW(std::out_of_range("Position of alternative allele is out of range of alleles in ALT column"));
    }

    if (!isInDNA5(record.alt[altPos])) {
        SEQAN_THROW(std::invalid_argument("Alternative allele is not a valid DNA base"));
    }

    seq[record.beginPos] = record.alt[altPos];
}

size_t getEndPos(const size_t beginPos, const std::unordered_map<std::string, std::string> & infoMap) {
    size_t endPos;

    auto it = infoMap.find("END");
    if (it == infoMap.end()) {
        it = infoMap.find("SVLEN");
        if (it == infoMap.end()) {
            throw ParseError("ERROR: Duplication record does not contain END or SVLEN in INFO field.");
        } else {
            endPos = beginPos + std::stoi(it->second)+1;
        }
    } else {
        endPos = std::stoi(it->second);

        if (endPos < beginPos) {
            throw std::out_of_range("ERROR: The value of key END in INFO column is smaller than the begin position.");
        }

    }

    return endPos;
}

size_t getTarget(CharString & targetContig, const std::unordered_map<std::string, std::string> & infoMap) {
    auto it = infoMap.find("TARGETPOS");

    if (it == infoMap.end()) {
        throw ParseError("ERROR: Duplication record does not contain TARGETPOS in INFO field.");
    }

    // split it.second by ':'
    std::istringstream iss{it->second};
    std::string contig, pos;
    std::getline(iss, contig, ':');
    std::getline(iss, pos);

    if (contig.empty() || pos.empty()) {
        throw ParseError("ERROR: Value of key TARGETPOS in INFO column is not of form contig:position.");
    }

    targetContig = contig;

    return std::stoi(pos) - 1;
}

template <typename TString>
bool getLine(std::istream & is, TString & line, const std::function<bool(char)> & delimiter, char & c) {
    while (is.get(c)) {
        if (delimiter(c)) {
            return true;
        }
        append(line, c);
    }
    return !empty(line);
}

template <typename TString>
bool getLine(std::istream & is, TString & line, const std::function<bool(char)> & delimiter) {
    char c{};
    return getLine(is, line, delimiter, c);
}

template <typename TString>
bool getLine(std::istream & is, TString & line) {
    char c{};
    while (is.get(c)) {
        append(line, c);
    }
    return !empty(line);
}

struct InDel;
struct Inv;
struct PosMap;
struct DupTarget;
struct Bnd;

struct ApplyVisitor {
    virtual ~ApplyVisitor() = default;
    virtual void visit(class InDel & inDel) = 0;
    virtual void visit(class Inv & inversion) = 0;
    virtual void visit(class PosMap & duplication) = 0;
    virtual void visit(class DupTarget & duplicationTarget) = 0;
    virtual void visit(class Bnd & breakEnd) = 0;
};

struct GetVisitor {
    virtual ~GetVisitor() = default;
    virtual void visit(const PosMap & duplication) = 0;
};

struct Var {
    virtual ~Var() = default;

    virtual void accept(ApplyVisitor& visitor) = 0;
    virtual void accept(GetVisitor& visitor) const = 0;
};

struct InDel : public Var {
    Dna5String ref{};
    Dna5String alt{};

    InDel(const VcfRecord & record, const unsigned altIdx) : ref(record.ref) {
        auto beginIt = begin(record.alt);

        for (size_t i = 1; i < altIdx; ++i) {
            beginIt = std::find(beginIt, end(record.alt), ',');
            if (beginIt == end(record.alt)) {
                throw ParseError("ERROR: Haplotype info was out of range.");
            }
            ++beginIt;
        }

        auto endIt = std::find(beginIt, end(record.alt), ',');

        append(alt, infix(record.alt, beginIt, endIt));
    }

    void accept(ApplyVisitor& visitor) override {
        visitor.visit(*this);
    }

    void accept(GetVisitor& visitor) const override {
    }

    friend std::ostream & operator<<(std::ostream & os, const InDel & inDel) {
        os << "\tInDel: ref=" << inDel.ref << ", alt=" << inDel.alt << "\n";
        return os;
    }

    [[nodiscard]] size_t apply(Dna5String & hapSeq) const {
        append(hapSeq, this->alt);

        std::cout << *this;

        return length(this->ref);
    }
};

struct Bnd : public Var {
    size_t targetPos{};
    CharString targetContig{};
    Dna5String alt{};
    bool isForward{};
    bool isRight{};

    explicit Bnd(const VcfRecord & record) {
        std::istringstream iss{toCString(record.alt)};
        std::function<bool(char)> delimiter = [](char c) { return c == '[' || c == ']'; };
        char c{};

        // I don't use strSplit because I need to know if the first character is [, ] or :
        if (!getLine(iss, alt, delimiter, c)) {
            throw ParseError("ERROR: Breakend record dose not contain [ or ].");
        }

        this->isForward = c == '[';
        this->isRight = empty(this->alt);

        if (!getLine(iss, this->targetContig, [](char c) { return c == ':'; })) {
            throw ParseError("ERROR: Breakend record does not contain : in ALT field.");
        }

        CharString pos{};

        if (!getLine(iss, pos, delimiter)) {
            throw ParseError("ERROR: Breakend record does not contain [ or ] in ALT field.");
        }

        this->targetPos = lexicalCast<size_t>(pos)-1;

        if (empty(this->alt) && !getLine(iss, this->alt)) {
            throw ParseError("ERROR: Breakend record does not contain any alt sequence after [. Must contain at least one base!");
        }
    }

    void accept(ApplyVisitor& visitor) override {
        visitor.visit(*this);
    }

    void accept(GetVisitor& visitor) const override {
    }

    friend std::ostream & operator<<(std::ostream & os, const Bnd & bnd) {
        os << "\tBnd: targetContig=" << bnd.targetContig << ", targetPos=" << bnd.targetPos << ", alt=" << bnd.alt << "\n";
        return os;
    }

    [[nodiscard]] size_t apply(Dna5String & hapSeq) const {
        append(hapSeq, this->alt);
        std::cout << * this;
        return targetPos;
    }
};

struct PosMap : public Var {
    size_t hapPos{};

    void accept(ApplyVisitor& visitor) override {
        visitor.visit(*this);
    }

    void accept(GetVisitor& visitor) const override {
        visitor.visit(*this);
    }

    friend std::ostream & operator<<(std::ostream & os, const PosMap & posMap) {
        os  << "\tPosMap: hapPos= " << posMap.hapPos << "\n";
        return os;
    }

    [[nodiscard]] size_t apply(const size_t currentLength) {
        this->hapPos = currentLength;
        std::cout << *this;
        return 0;
    }
};

struct Get : public GetVisitor {
    size_t hapPos{};

    void visit(const PosMap & posMap) override {
        hapPos = posMap.hapPos;
    }
};

struct Inv : public Var {
    size_t beginPos{};

    explicit Inv(const size_t beginPos) {
        this->beginPos = beginPos;
    }

    void accept(ApplyVisitor& visitor) override {
        visitor.visit(*this);
    }

    void accept(GetVisitor& visitor) const override {
    }

    friend std::ostream & operator<<(std::ostream & os, const Inv & inv) {
        os << "\tInv: beginPos= " << inv.beginPos << "\n";
        return os;
    }

    [[nodiscard]] size_t apply(Dna5String& hapSeq, std::map<size_t, std::unique_ptr<Var>>& vars) const {
        Get get{};

        const auto it = vars.find(beginPos);

        if (it == vars.end()) {
            throw std::out_of_range("ERROR: Could not find PosMap for inversion.");
        }

        it->second->accept(get);
        std::cout << length(infix(hapSeq, get.hapPos, length(hapSeq))) << "\n";
        std::cout << length(hapSeq) << "\n";
        reverseComplement(infix(hapSeq, get.hapPos, length(hapSeq)));

        std::cout << *this;

        return 0;
    }
};

struct DupTarget : public Var {
    size_t sourceBeginPos{};
    size_t sourceEndPos{};

    CharString targetContig{};

    DupTarget(size_t sourcePos, size_t endPos) :
        sourceBeginPos(sourcePos), sourceEndPos(endPos) {
    }

    void accept(ApplyVisitor& visitor) override {
        visitor.visit(*this);
    }

    void accept(GetVisitor& visitor) const override {
    }

    friend std::ostream & operator<<(std::ostream & os, const DupTarget & dupTarget) {
        os << "\tDupTarget: sourceBeginPos=" << dupTarget.sourceBeginPos << "sourceEndPos=" << dupTarget.sourceEndPos << ", targetContig=" << dupTarget.targetContig << "\n";
        return os;
    }

    [[nodiscard]] size_t apply(Dna5String& hapSeq, std::map<size_t, std::unique_ptr<Var>>& vars) const {

        Get getBegin{};
        Get getEnd{};

        const auto itBegin = vars.find(sourceBeginPos);
        const auto itEnd = vars.find(sourceEndPos);

        if (itBegin == vars.end() || itEnd == vars.end()) {
            throw std::out_of_range("ERROR: Could not find PosMap for duplication source begin.");
        }

        itBegin->second->accept(getBegin);
        itEnd->second->accept(getEnd);

        append(hapSeq, infix(hapSeq, getBegin.hapPos, getEnd.hapPos));

        std::cout << *this;
        return 0;
    }
};

template <typename It>
struct Apply : public ApplyVisitor {
    Dna5String& hapSeq;
    Dna5String& seq;
    std::map<size_t, std::unique_ptr<Var>>& vars;
    It& i;

    size_t pos{};
    bool visitedBnd{};

    explicit Apply(Dna5String& hapSeq, Dna5String& seq, std::map<size_t, std::unique_ptr<Var>>& vars, It& i) :
                    hapSeq(hapSeq), seq(seq), vars(vars), i{i} {}

    void visit(InDel& inDel) override {
        this->pos = inDel.apply(this->hapSeq);
    }

    void visit(Inv& inversion) override {
        pos = inversion.apply(this->hapSeq, this->vars);
    }

    void visit(PosMap& duplicationSource) override {
        pos = duplicationSource.apply(length(hapSeq));
    }

    void visit(DupTarget& duplicationTarget) override {
        pos = duplicationTarget.apply(this->hapSeq,this->vars);
    }

    void visit(Bnd& breakEnd) override {
        pos = breakEnd.apply(this->hapSeq);
        this->visitedBnd = true;
    }
};

bool loadSV(std::map<size_t, std::unique_ptr<Var>> & vars,
            std::unordered_map<std::string, std::string> & infoMap, const VcfRecord & record) {

    if (infoMap["SVTYPE"] == "INV") {
        vars.emplace(record.beginPos+1, std::make_unique<PosMap>());
        vars.emplace(getEndPos(record.beginPos, infoMap), std::make_unique<Inv>(record.beginPos+1));
    } else if (infoMap["SVTYPE"] == "DUP") {
        ////TODO: handle tandem repeats without target pos
        ////TODO: check if target is out of range

        CharString targetContig;
        size_t endPos = getEndPos(record.beginPos, infoMap);

        //// TODO: why dose inversion and duplications point to the last ref and not the first mutated base?
        vars.emplace(record.beginPos+1, std::make_unique<PosMap>());
        vars.emplace(endPos, std::make_unique<PosMap>());
        vars.emplace(getTarget(targetContig, infoMap), std::make_unique<DupTarget>(record.beginPos+1, endPos));
    } else if (infoMap["SVTYPE"] == "BND") {
        ////TODO: check if mate is out of range
        vars.emplace(record.beginPos, std::make_unique<Bnd>(record));
    }
    else {
        std::cerr << "ERROR: This kind of variant is not supported. " << record.id << "\n";
        return false;
    }
    return true;
}

void applyVars(Dna5String & hapSeq, Dna5String & seq, std::map<size_t, std::unique_ptr<Var>> & vars) {
    size_t beginPos{0};
    size_t endPos;
    auto it = vars.begin();

    Apply visitor{hapSeq, seq, vars, it};

    for (; it != vars.end(); ) {
        endPos = it->first;
        // Copy the sequence between the last variant and the current one
        //// TODO: handle overlapping variants (for example when a duplication source is overlapping with a deletion)
        append(hapSeq, infix(seq, beginPos, endPos));
        // Apply the variance
        it->second->accept(visitor);
        // Load new begin Pos

        // find next variant
        if (!visitor.visitedBnd) {
            beginPos = endPos + visitor.pos;
            ++it;
        } else {
            beginPos = visitor.pos;
            it = vars.find(beginPos);

            if (it->first != beginPos) {
                throw ParseError("ERROR: Left side of breakend dose not match with right side!\n");
            } else {
                ++it;
            }
            visitor.visitedBnd = false;
        }
    }

    endPos = length(seq);
    append(hapSeq, infix(seq, beginPos, endPos));
}

void loadRecord(VcfFileIn & vcfIn, const unsigned haploIdx, const unsigned contigIdx,
                Dna5String & seq, std::map<size_t, std::unique_ptr<Var>> & vars) {
    VcfRecord record{};

    readRecord(record, vcfIn);

    if (record.rID < contigIdx) {
        throw ParseError("VCF file is not sorted by position.");
    } else if (record.rID > contigIdx) {
        // handle contig change
    } else {
        if (unsigned altIdx = getAltIdx(record, haploIdx); altIdx != 0) {
            std::unordered_map<std::string, std::string> infoMap{};
            getInfo(infoMap, record);

            if (infoMap.find("SVTYPE") == infoMap.end()) {
                if (isSNP(record)) {
                    applySNP(seq, record, altIdx);
                } else if (isInDel(record)) {
                    vars.emplace(record.beginPos, std::make_unique<InDel>(record, altIdx));
                } else {
                    //// TODO: is there a better way to print out CharStrings?
                    throw ParseError("ERROR: This kind of variant is not supported. " + std::string(toCString(record.id)));
                }
            } else {
                loadSV(vars, infoMap, record);
            }
        }
    }
}

bool createHaplotypeSequence(SeqFileOut & out, const FaiIndex & ref, VcfFileIn & vcfIn, const unsigned haploIdx) {

    bool res = true;

    for (const CharString & contigName : ref.seqNameStore) {

        unsigned contigIdx = 0;
        getIdByName(contigIdx, ref, contigName);

        std::map<size_t, std::unique_ptr<Var>>  vars{};

        Dna5String seq{};
        readSequence(seq, ref, contigIdx);

        while (!atEnd(vcfIn)) {
            try {
                loadRecord(vcfIn, haploIdx, contigIdx, seq, vars);
            } catch (ParseError const &e) {
                std::cerr << "ERROR: Could not parse record. " << e.what() << "\n";
                res = false;
                break;
            } catch (IOError const &e) {
                std::cerr << "ERROR: Could not read record. " << e.what() << "\n";
                res = false;
                break;
            }
        }

        try {
            Dna5String hapSeq{};
            applyVars(hapSeq, seq, vars);
            ////TODO: write out directly to file
            writeRecord(out, contigName, hapSeq);
        } catch (IOError const &e) {
            std::cerr << "ERROR: Could not write haplotype sequence to file. " << e.what() << "\n";
            res = false;
        }
    }

    return res;
}
