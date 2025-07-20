#include <vector>
#include <array>
#include <iostream>
#include <stdexcept>
#include <cmath>


#define nClrs 8
#define nSlots 4

class clr_t{
private:
    unsigned int clr;
public:
    clr_t(unsigned int clr = 0)
        :clr(clr) {
            if (nClrs < clr) clr = -1;
        }

    operator unsigned int() const {
        return clr;
    }

    clr_t& operator++(){
        ++clr;
        return *this;
    }

    std::string toString() const{
        if (clr >= nClrs) return "x";
        return std::string(1, (char)('A' + clr));
    }
};

std::ostream& operator<<(std::ostream& out, const clr_t& c)
{
   return out << c.toString();
}
 


typedef std::vector<clr_t> seq_t;

std::string seqToString(seq_t prefix){
    std::string res("[");
    size_t length(prefix.size());
    for (size_t i(0); i < nSlots; ++i){
        if (i < length) res += prefix[i].toString();
        else res += "*";
    }
    res += "]";
    return res;
}

void myPrint(seq_t prefix){
    size_t length(prefix.size());
    std::cout << "[";
    for (size_t i(0); i < nSlots; ++i){
        if (i < length) std::cout << prefix[i];
        else std::cout << "*";
    }
    std::cout << "]";
}

class res_t{
private:
    unsigned int correct;
    unsigned int misplaced;
public:
    res_t(unsigned int correct = 0, unsigned int misplaced = 0)
        :correct(correct), misplaced(misplaced) {
            if (nSlots < correct) throw std::invalid_argument("Error: res_t out of range");
            if (nSlots < misplaced) throw std::invalid_argument("Error: res_t out of range");
    }
    
    unsigned int getCorrect() const {
        return correct;
    }

    unsigned int getMisplaced() const {
        return misplaced;
    }

    bool operator==(const res_t& other) const {
        return (correct == other.correct && misplaced == other.misplaced);
    }
};

class resCoun_t;

class rule_t{
private:
    seq_t seq;
    res_t res;
    bool hasRes;
public:
    rule_t(seq_t const& seq, res_t const& res)
        :seq(seq), res(res), hasRes(true){
            if (nSlots != seq.size()) throw std::invalid_argument("Error: rule_t seq size mismatch");
    }

    rule_t(seq_t const& seq)
        :seq(seq), hasRes(false){
            if (nSlots != seq.size()) throw std::invalid_argument("Error: rule_t seq size mismatch");
    }

    rule_t(seq_t const& seq, seq_t const& secretCode)
        :seq(seq), hasRes(true){
            if (nSlots != seq.size()) throw std::invalid_argument("Error: rule_t seq size mismatch");
            if (nSlots != secretCode.size()) throw std::invalid_argument("Error: rule_t secretCode size mismatch");
            res = getResFrom(secretCode);
        }

    bool doesAccept(seq_t const& prefix) const {
        if (prefix.size() != seq.size()) throw std::invalid_argument("Error: rule_t seq size mismatch");
        if (not hasRes) throw std::logic_error("Error: rule_t has not result but doesAccept is called");
        res_t otherRes(getResFrom(prefix));
        return (otherRes == res);
    }

    res_t getResFrom(seq_t prefix) const {
        // std::cout << std::endl;
        // std::cout << "START" << std::endl;
        // std::cout << "Calculating result for prefix: " << seqToString(prefix) << std::endl;
        // std::cout << "                               " << seqToString(seq) << std::endl;
        
        unsigned int correct(0);
        std::array<bool, nSlots> used({});
        for (size_t i(0); i < nSlots; ++i){
            if (prefix[i] == seq[i]) {
                ++correct;
                used[i] = true;
                prefix[i] = -1; // Mark as used
            }
        }

        // std::cout << "Calculating result for prefix: " << seqToString(prefix) << std::endl;
        // std::cout << "                               " << seqToString(seq) << std::endl;

        unsigned int misplaced(0);
        for (size_t i(0); i < nSlots; ++i){
            if (not used[i]){
                for (size_t j(0); j < nSlots; ++j){
                    if (i != j){
                        // std::cout << "Comparing " << i << ", " << j << " --> " << seq[i] << " with " << prefix[j] << std::endl;
                        if (prefix[j] == seq[i]) {
                            ++misplaced;
                            prefix[j] = -1;
                            j = nSlots; // Break inner loop
                            // std::cout << "FOUND" << std::endl;
                            // std::cout << "Calculating result for prefix: " << seqToString(prefix) << std::endl;
                            // std::cout << "                               " << seqToString(seq) << std::endl;
                        }
                    }
                }
            }
        }
        // std::cout << "--> " << correct << " correct, " << misplaced << " misplaced" << std::endl;
        return res_t(correct, misplaced);
    }

    std::string toString() const{
        if (hasRes) return seqToString(seq) + "(" + std::to_string(res.getCorrect()) + ", " + std::to_string(res.getMisplaced()) + ")"; 
        return seqToString(seq) + "(?, ?)"; 
    }
};

std::ostream& operator<<(std::ostream& out, const rule_t& rule)
{
   return out << rule.toString();
}

class resCoun_t{
private:
    std::array<std::array<unsigned int, nSlots+1>, nSlots+1> counts;
public:
    resCoun_t() 
        :counts() {}

    resCoun_t(res_t res)
        :counts() {
            if (nSlots < res.getCorrect()) throw std::invalid_argument("Error: resCoun_t out of range");
            if (nSlots < res.getMisplaced()) throw std::invalid_argument("Error: resCoun_t out of range");
            counts[res.getCorrect()][res.getMisplaced()] = 1;
    }

    resCoun_t(unsigned int correct, unsigned int misplaced)
        :counts() {
            if (nSlots < correct) throw std::invalid_argument("Error: resCoun_t out of range");
            if (nSlots < misplaced) throw std::invalid_argument("Error: resCoun_t out of range");
            counts[correct][misplaced] = 1;
    }

    resCoun_t operator+=(const resCoun_t& other) {
        for (unsigned int i(0); i <= nSlots; ++i){
            for (unsigned int j(0); j <= nSlots; ++j){
                counts[i][j] += other.counts[i][j];
            }
        }
        return *this;
    }

    std::string toString() const {
        std::string res = "Result counts:\n";
        for (unsigned int i(0); i <= nSlots; ++i){
            for (unsigned int j(0); j <= nSlots; ++j){
                if (0 < counts[i][j]) {
                    res += "Correct: " + std::to_string(i) + ", Misplaced: " + std::to_string(j) + " -> Count: " + std::to_string(counts[i][j]) + "\n";
                }
            }
        }
        return res;
    }
};


std::ostream& operator<<(std::ostream& out, const resCoun_t& res)
{
   return out << res.toString();
}

class explorationResul_t {
    unsigned int nPoss;
public:
    explorationResul_t(unsigned int nPoss = 0)
        :nPoss(nPoss) {}

    explorationResul_t operator+=(const explorationResul_t& other) {
        nPoss += other.nPoss;
        return *this;
    }

    std::string toString() const {
        return "Number of possibilities: " + std::to_string(nPoss);
    }

    unsigned int getCount() const {
        return nPoss;
    }
};

std::ostream& operator<<(std::ostream& out, const explorationResul_t& res)
{
   return out << res.toString();
}



class Node;

class Node{
protected:
    bool isDead;
public:
    Node()
        :isDead(0) {}
    virtual Node* getChildren(clr_t clr) const = 0;
    virtual explorationResul_t exploreAllPoss(seq_t prefix, bool printAll) const = 0;
    virtual bool filterFor(seq_t prefix, rule_t const& rule) = 0;
    virtual resCoun_t statsFor(seq_t prefix, rule_t const& guess) const = 0;
};


class LeafNode : public Node{
public:
    LeafNode()
    :Node() {}

    virtual Node* getChildren(clr_t) const override final{
        throw std::logic_error("Error - LeafNode::getChildren called");
    }

    virtual explorationResul_t exploreAllPoss(seq_t prefix, bool printAll) const override final{
        if (isDead) return explorationResul_t(0);
        if (printAll){
            std::cout << "* ";
            myPrint(prefix);
            std::cout << std::endl;
        }
        return explorationResul_t(1);
    }

    virtual bool filterFor(seq_t prefix, rule_t const& rule) override final{
        if (isDead) return false;
        isDead = not rule.doesAccept(prefix);
        return not isDead;
    }

    virtual resCoun_t statsFor(seq_t prefix, rule_t const& guess) const override final{
        if (isDead) return resCoun_t();

        res_t res(guess.getResFrom(prefix));
        resCoun_t resCount(res);
        return resCount;
    }
};


class ParentNode : public Node{
private:
    std::array<Node*, nClrs> children; // Should be vector for post-compilation input by user
public:
    ParentNode(unsigned int depth)
        :Node() {
            for (clr_t clr(0); clr < nClrs; ++clr){
                if (depth +1 < nSlots){ //WARNING
                    children[clr] = new ParentNode(depth+1);
                } else {
                    children[clr] = new LeafNode();
                }
            }
        }

    virtual Node* getChildren(clr_t clr) const override final{
        return children[clr];
    }

    virtual explorationResul_t exploreAllPoss(seq_t prefix, bool printAll) const override final{
        explorationResul_t resExp(0);
        if (isDead) return resExp;
        seq_t newPrefix(prefix);
        newPrefix.push_back(0);

        for (clr_t clr(0); clr < nClrs; ++clr){
            newPrefix.back() = clr;
            resExp += getChildren(clr)->exploreAllPoss(newPrefix, printAll);
        }
        return resExp;
    }

    virtual bool filterFor(seq_t prefix, rule_t const& rule) override final{
        if (isDead) return false;

        seq_t newPrefix(prefix);
        newPrefix.push_back(0);

        bool hasLivingChildren(false);
        for (clr_t clr(0); clr < nClrs; ++clr){
            newPrefix.back() = clr;
            bool isAlive(getChildren(clr)->filterFor(newPrefix, rule));
            hasLivingChildren = hasLivingChildren or isAlive;
        }
        if (!hasLivingChildren) isDead = true;
        return hasLivingChildren;
    }

    virtual resCoun_t statsFor(seq_t prefix, rule_t const& guess) const override final{
        resCoun_t resCount;
        if (isDead) return resCount;

        seq_t newPrefix(prefix);
        newPrefix.push_back(0);

        for (clr_t clr(0); clr < nClrs; ++clr){
            newPrefix.back() = clr;
            resCoun_t thisRes(getChildren(clr)->statsFor(newPrefix, guess));
            resCount += thisRes;
        }
        return resCount;
    }
};


class RootNode : public ParentNode{
public:
    RootNode()
    :ParentNode(0) {}

    explorationResul_t exploreAllPoss(bool printAll = true) const{
        if (printAll) std::cout << "Printing all possibilities:" << std::endl;
        seq_t emptySequence({});
        explorationResul_t res(ParentNode::exploreAllPoss(emptySequence, printAll));
        if (printAll) std::cout << res << std::endl;
        return res;
    }

    void filterFor(rule_t const& rule){
        seq_t emptySequence({});
        std::cout << "Filtering for information " << rule << std::endl;
        std::cout << std::endl;
        ParentNode::filterFor(emptySequence, rule);
    }

    void statsFor(seq_t const& guess) const{
        seq_t emptySequence({});
        rule_t semiRule(guess);
        resCoun_t resCount(ParentNode::statsFor(emptySequence, semiRule));
        std::cout << "Analysis for guess ";
        myPrint(guess);
        std::cout << std::endl;
        std::cout << resCount << std::endl;
    }
};


void generateAllStartingPatterns(seq_t prefix, size_t maxRepetitions, std::vector<seq_t>& allPatterns){
    // std::cout << "Generating all starting patterns for prefix: " << seqToString(prefix) << " - " << maxRepetitions << std::endl;
    size_t length(prefix.size());
    clr_t firstAllowed((length == 0 ? 0 : prefix.back()+1));
    if (nClrs <= firstAllowed) return;
    for (size_t n(1); n <= maxRepetitions && length + n <= nSlots; ++n){
        seq_t newPrefix(prefix);
        for (size_t i(0); i < n; ++i){
            newPrefix.push_back(firstAllowed);
        }
        if (length + n == nSlots) {
            allPatterns.push_back(newPrefix);
        } else {
            generateAllStartingPatterns(newPrefix, std::min(maxRepetitions, n), allPatterns);
        }
    }
    // std::cout << "END       all starting patterns for prefix: " << seqToString(prefix) << " - " << maxRepetitions << std::endl;
}


void exploreFirstMove(RootNode const& possibilities){
    std::cout << "Let's explore the first move together." << std::endl;
    if (possibilities.exploreAllPoss(false).getCount() != std::pow(nClrs, nSlots)){
        throw std::logic_error("Error: First move exploration does not yield all possibilities");
    }
    std::cout << "There are " << std::pow(nClrs, nSlots) << " possibilities." << std::endl;

    std::vector<seq_t> allPatterns({});
    generateAllStartingPatterns(seq_t({}), nSlots, allPatterns);

    for (size_t i(0); i < allPatterns.size(); ++i){
        seq_t const& pattern(allPatterns[i]);
        std::cout << "Exploring pattern " << i+1 << "/" << allPatterns.size() << ": " << seqToString(pattern) << std::endl;
        possibilities.statsFor(pattern);
    }
}



int main(){
    RootNode possibilities;

    exploreFirstMove(possibilities);
    /*
    possibilities.exploreAllPoss();
    std::cout << std::endl;

    seq_t guess_allZero({0, 0, 0, 0});

    possibilities.statsFor(guess_allZero);
    std::cout << std::endl;

    rule_t rule_TwoZero(guess_allZero, res_t(2, 0));
    possibilities.filterFor(rule_TwoZero);
    std::cout << std::endl;

    possibilities.exploreAllPoss();
    std::cout << std::endl;


    seq_t guess_allOnes({1, 1, 1, 1});
    possibilities.statsFor(guess_allOnes);
    std::cout << std::endl;

    rule_t rule_OneOne(guess_allOnes, res_t(1, 0));
    possibilities.filterFor(rule_OneOne);
    std::cout << std::endl;


    possibilities.exploreAllPoss();
    std::cout << std::endl;


    seq_t guess_AABB({0, 0, 1, 1});
    possibilities.statsFor(guess_AABB);
    std::cout << std::endl;

    rule_t rule_AABB(guess_AABB, res_t(2, 1));
    possibilities.filterFor(rule_AABB);
    std::cout << std::endl;


    possibilities.exploreAllPoss();
    std::cout << std::endl;


    seq_t guess_ACAB({0, 2, 0, 1});
    possibilities.statsFor(guess_ACAB);
    std::cout << std::endl;

    rule_t rule_ACAB(guess_ACAB, res_t(2, 2));
    possibilities.filterFor(rule_ACAB);
    std::cout << std::endl;

    possibilities.exploreAllPoss();
    std::cout << std::endl;
    */


    return 0;
}