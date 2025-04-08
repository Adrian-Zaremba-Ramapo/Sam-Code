import copy

var('a'); var('b'); var('n')

def listMultiply(numList1, numList2):
    result = []
    for num1 in numList1:
        for num2 in numList2:
            result.append(num1 * num2)
    return result

def getAllFactorsRecursive(factorList):
    if len(factorList) == 0:
        return [1]
    firstElement = factorList.pop(0)
    firstFactor = firstElement[0]
    multiplicity = firstElement[1]
    firstFactorList = []
    for i in range(multiplicity + 1):
        firstFactorList.append(firstFactor^i)
    return listMultiply(firstFactorList, getAllFactorsRecursive(factorList))

def getAllFactors(num):
    return sorted(getAllFactorsRecursive(list(factor(num))))

def getPossibleGCDs(index, value):
    if SR(value).degree(n) == 0:
        return [1]
    indexNCoefficient = ZZ(index.coefficient(n, 1).coefficient(a, 1))
    valueNCoefficient = ZZ(value.coefficient(n, 1).coefficient(a, 1))
    coefficientLCM = ZZ(lcm(indexNCoefficient, valueNCoefficient))
    indexMultiplier = coefficientLCM // indexNCoefficient
    valueMultiplier = coefficientLCM // valueNCoefficient
    difference = index * indexMultiplier - value * valueMultiplier
    if difference != 0:
        return getAllFactors(ZZ(difference))
    coefficientGCD = ZZ(gcd(indexNCoefficient, valueNCoefficient))
    return [index / (indexNCoefficient // coefficientGCD)]

def recursiveSearch(index, value, depth, currPath, allPaths, prevAdditionCount, prevDivisionCount):
    if str(value) == "15*a*n + 13":
        max(1, 1)
    currPath.append((index, value))
    if depth <= 1:
        allPaths.append(currPath)
        return
    nextIndex = index + 1
    indexACoefficient = ZZ(SR(nextIndex).coefficient(n, 1).coefficient(a, 1))
    indexBCoefficient = ZZ(SR(nextIndex).coefficient(n, 0).coefficient(b, 1))
    indexConstant = ZZ(SR(nextIndex).coefficient(n, 0).coefficient(b, 0))
    indexGCD = gcd([indexACoefficient, indexBCoefficient, indexConstant])
    valueACoefficient = ZZ(SR(value).coefficient(n, 1).coefficient(a, 1))
    valueBCoefficient = ZZ(SR(value).coefficient(n, 0).coefficient(b, 1))
    valueConstant = ZZ(SR(value).coefficient(n, 0).coefficient(b, 0))
    valueGCD = gcd([valueACoefficient, valueBCoefficient, valueConstant])
    totalGCD = gcd(indexGCD, valueGCD)
    possibleGCDs = getPossibleGCDs(nextIndex, value)
    for possibleGCD in possibleGCDs:
        if possibleGCD not in ZZ:
            nextValue = nextIndex / possibleGCD
            nextCurrPath = copy.copy(currPath)
            recursiveSearch(nextIndex, nextValue, depth - 1, nextCurrPath, allPaths, 0, 1)
            continue
        if possibleGCD % totalGCD != 0:
            continue
        if possibleGCD == 1:
            if prevAdditionCount < 3:
                nextValue = value + nextIndex + 2
                nextCurrPath = copy.copy(currPath)
                recursiveSearch(nextIndex, nextValue, depth - 1, nextCurrPath, allPaths, prevAdditionCount + 1, 0)
        else:
            if prevDivisionCount == 0:
                divisionFactor, inverse, _ = xgcd(indexBCoefficient, possibleGCD)
                if indexConstant % divisionFactor == 0:
                    newIndexBCoefficient = possibleGCD // divisionFactor
                    newIndexConstant = indexConstant // divisionFactor
                    newIndexBTerm = newIndexBCoefficient * b - (newIndexConstant * inverse % possibleGCD)
                    newIndexACoefficient = ZZ(SR(nextIndex).coefficient(n, 1).coefficient(a, 1))
                    newIndexATerm = (lcm(newIndexACoefficient, possibleGCD) // newIndexACoefficient) * a
                    updatedNextIndex = nextIndex(b = newBTerm)(a = newATerm)
                    
                    divisionFactor, inverse, _ = xgcd(newIndexBCoefficient * valueBCoefficient, possibleGCD)
                    if valueConstant % divisionFactor == 0:
                        newValueBCoefficient = possibleGCD // divisionFactor
                        newValueConstant = indexConstant // divisionFactor
                        newValueBTerm = newIndexBCoefficient * b - (newIndexConstant * inverse % possibleGCD)
                        newValueACoefficient = ZZ(SR(nextIndex).coefficient(n, 1).coefficient(a, 1))
                        newIndexATerm = (lcm(newIndexACoefficient, possibleGCD) // newIndexACoefficient) * a

                        
                        nextValue = updatedNextIndex / possibleGCD
                        nextCurrPath = []
                        for termPair in currPath:
                            nextCurrPathIndex = SR(termPair[0])(b = newIndexBTerm)(a = newATerm)
                            nextCurrPathValue = SR(termPair[1])(b = newIndexBTerm)(a = newATerm)
                            nextCurrPath.append((nextCurrPathIndex, nextCurrPathValue))
                        recursiveSearch(updatedNextIndex, nextValue, depth - 1, nextCurrPath, allPaths, 0, 1)



def centerPad(string, n):
    padLength = n - len(string)
    leftPadLength = padLength // 2
    rightPadLength = padLength - leftPadLength
    return ' ' * leftPadLength + string + ' ' * rightPadLength

def printNElements(n, start, table, termLength):
    for i in range(start, min(start + n, len(table))):
        element = table[i]
        print(centerPad(str(element[0]), termLength), end = '')
        print(' | ', end = '')
    print()
    for i in range(start, min(start + n, len(table))):
        element = table[i]
        print(centerPad(str(element[1]), termLength), end = '')
        print(' | ', end = '')
    print()

def printTable(table, termLength, termsPerLine):
    for i in range(len(table) // termsPerLine):
        printNElements(termsPerLine, i * termsPerLine, table, termLength)
        print()
    if len(table) % termsPerLine != 0:
        printNElements(len(table) % termsPerLine, len(table) - len(table) % termsPerLine, table, termLength)

allPaths = []
depth = 11
recursiveSearch(12*a*n, 2, depth, [], allPaths, 0, 0)
print('*' * 100)
for path in allPaths:
    printTable(path, 18, 5)
    print('*' * 100)

                    
