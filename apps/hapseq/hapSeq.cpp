//
// Created by dominik on 20.06.24.
//

#include "variance.cpp"

using namespace seqan2;

int main (int argc, char *argv[]) {
    ConsensusOptions options{};

    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != ArgumentParser::PARSE_OK) {
        return res == ArgumentParser::PARSE_ERROR;
    }

    FaiIndex ref{};

    if (!loadOrBuildFAI(ref, options.referencePath)) {
        std::cerr << "Could not load or build FAI index for " << options.referencePath << "\n";
        return 1;
    }

    VcfFileIn vcfIn{}; VcfHeader vcfHeader{};

    if (!openAndLoadVCFHeader(vcfHeader, vcfIn, options.variancePath)) {
        std::cerr << "Could not open and load VCF header for " << options.variancePath << "\n";
        return 1;
    }

    const StringSet<CharString> & contigNameSet{contigNames(context(vcfIn))};
    const StringSet<CharString> & sampleNameSet{sampleNames(context(vcfIn))};

    for (size_t i = 0; i < length(contigNameSet); ++i) {
        if (ref.seqNameStore[i] != contigNameSet[i]) {
            std::cerr << "Reference and variance file do not match! The " << i << " contig is different\n";
            std::cerr << "Reference: " << ref.seqNameStore[i] << " Variance: " << contigNameSet[i] << "\n";
            return 1;
        }
    }

    for (const CharString & haplotype : options.haplotypeNames) {
        std::filesystem::path outPath =
                options.outdirPath / (static_cast<std::string>(toCString(haplotype)) + ".fa");
        SeqFileOut out{};

        unsigned haplotypeIdx{};
        getIdByName(haplotypeIdx, sampleNameSet,  haplotype);

        if (!open(out, outPath.c_str())) {
            std::cerr << "Could not open output file " << outPath << std::endl;
            return 1;
        }

        if (!createHaplotypeSequence(out, ref, vcfIn, haplotypeIdx)) {
            std::cerr << "Could not create haplotype sequence for " << haplotype << "\n";
            return 1;
        }
    }

    return 0;
}