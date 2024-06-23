//
// Created by dominik on 20.06.24.
//

#include <seqan/vcf_io.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <filesystem>
#include <iostream>

using namespace seqan2;

struct ConsensusOptions {
    std::filesystem::path variancePath{};
    std::filesystem::path referencePath{};
    String<CharString> haplotypeNames{};
    std::filesystem::path outdirPath{};
};

void setupParser(ArgumentParser & parser) {
    setShortDescription(parser, "Generate consensus sequences from a reference and a variance file.");
    setVersion(parser, "0.0.1");
    setDate(parser, "June 2024");

    addUsageLine(parser,
                 R"([\fIOPTIONS\fP] "\fITEXT\fP")");
    addDescription(parser,
                   "This program creates a consensus sequences from a variance file in vcf format and reference"
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
                                             "Name of haplotype to create consensus sequence for. "
                                             "Can be given multiple times.",
                                             ArgParseArgument::STRING,
                                             "STR",
                                             true));

    addSection(parser, "Mandatory Output");

    addOption(parser, ArgParseOption("o",
                                             "outdir",
                                             "Path to directory, where the consensus files get saved "
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

    ArgumentParser parser("consensus");

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

bool createHaplotypeSequence(StringSet<Dna5String> & seqs, const FaiIndex & ref, VcfFileIn & vcfIn) {

    bool res = true;

    VcfRecord record{};

    while (!atEnd(vcfIn)) {
        try {
            readRecord(record, vcfIn);
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

    return res;
}

int main (int argc, char *argv[]) {

    ConsensusOptions options{};

    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != ArgumentParser::PARSE_OK) {
        return res == ArgumentParser::PARSE_ERROR;
    }

    FaiIndex ref{};

    if (!loadOrBuildFAI(ref, options.referencePath)) {
        return 1;
    }

    VcfFileIn vcfIn{}; VcfHeader header;

    if (!openAndLoadVCFHeader(header, vcfIn, options.variancePath)) {
        return 1;
    }

    const StringSet<CharString> & contigNameSet{contigNames(context(vcfIn))};

    for (size_t i = 0; i < length(contigNameSet); ++i) {
        if (ref.seqNameStore[i] != contigNameSet[i]) {
            std::cerr << "Reference and variance file do not match! The " << i << "contig is different\n";
            std::cerr << "Reference: " << ref.seqNameStore[i] << " Variance: " << contigNameSet[i] << "\n";
            return 1;
        }
    }

    for (const CharString & haplotype : options.haplotypeNames) {
        std::filesystem::path outPath =
                options.outdirPath / (static_cast<std::string>(toCString(haplotype)) + ".fa");
        SeqFileOut out{};

        if (!open(out, outPath.c_str())) {
            std::cerr << "Could not open output file " << outPath << std::endl;
            return 1;
        }

        StringSet<Dna5String> seqs{};
        resize(seqs, length(contigNameSet));

        for (size_t i = 0; i < length(contigNameSet); ++i) {
            reserve(seqs[i], sequenceLength(ref, i));
        }

        if (createHaplotypeSequence(seqs, ref, vcfIn)) {
            std::cerr << "Could not create haplotype sequence for " << haplotype << "\n";
            return 1;
        }

        try {
            writeRecords(out, contigNameSet, seqs);
        } catch (IOError const & e) {
            std::cerr << "ERROR: Could not write haplotype sequence to file. " << e.what() << "\n";
            return 1;
        }
    }

    return 0;
}