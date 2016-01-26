


class Amplicon:
    '''
    A complete Amplicon for use in EATING; the class contains methods for designing
    amplifying primers and describing the guides within it.

    properties:
    .start, .end: count of runs along chromosome
    .fiveprimeabut, .threeprimeabut: count of guides along chromosome that abut good guides
    .guides: list of good guides in run
    .length: permissible region distance
    .guidecount: number of guides in list (=len(self.guides))
    .genomename: genomename; must match BLAST database name for screening
    .chromosome
    .left_outside: the id of the guide to the left of the outside
    .left_inside: ...
    .permissible_start: absolute numbering vs chromosome
    .required_start_absolute: absolute numbering vs chromosome
    .required_start_relative
    .permissible_region
    '''
    def __init__(self, run, guidelist, genomename, chromosome):
        self.start = run[0]
        self.end = run[1]
        try:
            self.fiveprimeabut = guidelist[self.start-1]
            self.threeprimeabut= guidelist[self.end]
        except:
            self.fiveprimeabut = "0"
            self.threeprimeabut = "None"
        self.guides = guidelist[self.start:self.end]
        self.length = int(self.threeprimeabut.name)-int(self.fiveprimeabut.name)
        self.guidecount = self.end - self.start
        self.genomename = genomename
        self.chromosome = chromosome

        # Figure out how far from the desired amplicon primers can be designed
        # without bumping into the next (non-specific) primer.
        self.left_outside = self.fiveprimeabut.id[-1]
        self.left_inside = self.guides[0].id[-1]
        if self.left_outside == "F" and self.left_inside == "R":
            self.permissible_start = int(self.fiveprimeabut.name) + 10
            self.required_start_absolute = int(self.guides[0].name) +14
        elif self.left_outside == "R" and self.left_inside == "R":
            self.permissible_start = int(self.fiveprimeabut.name) + 1
            self.required_start_absolute = int(self.guides[0].name) +14
        elif self.left_outside == "R" and self.left_inside == "F":
            self.permissible_start = int(self.fiveprimeabut.name) + 1
            self.required_start_absolute = int(self.guides[0].name) +18
        elif self.left_outside == "F" and self.left_inside == "F":
            self.permissible_start = int(self.fiveprimeabut.name) + 10
            self.required_start_absolute = int(self.guides[0].name) +18
        else:
            print("error on left")
        # (fiveprimeabuttingguide, threeprimeabuttingguide), end-start, scores_and_details[start:end]))
        #self.right_inside = item[2][-1][1].id[-1]
        self.right_inside = self.guides[-1].id[-1]
        self.right_outside = self.threeprimeabut.id[-1]
        if self.right_outside == "F" and self.right_inside == "R":
            self.permissible_end = int(self.threeprimeabut.name) + 19
            self.required_end_absolute = int(self.guides[-1].name) + 2
        elif self.right_outside == "R" and self.right_inside == "F":
            self.permissible_end = int(self.threeprimeabut.name) + 10
            self.required_end_absolute = int(self.guides[-1].name) + 8
        elif self.right_outside == "R" and self.right_inside == "R":
            self.permissible_end = int(self.threeprimeabut.name) + 10
            self.required_end_absolute = int(self.guides[-1].name) + 2
        elif self.right_outside == "F" and self.right_inside == "F":
            self.permissible_end = int(self.threeprimeabut.name) + 19
            self.required_end_absolute = int(self.guides[-1].name) + 8
        else:
            print("error on right")

        self.permissible_region = self.chromosome[self.permissible_start:self.permissible_end]
        # Bounds that need to be included in PCR product :
        self.required_start_relative = self.required_start_absolute-self.permissible_start
        self.required_end_relative = self.required_end_absolute - self.permissible_start
        # Set up some other convenience stuff:
        self.permissible_region.name =str(self.fiveprimeabut.name)
        self.permissible_region.id =str(self.fiveprimeabut.name)
        self.permissible_region.description=str(self.guidecount)
        self.permissible_region.seq.alphabet = IUPACAmbiguousDNA()

    def __str__(self):
        return \
        "Maximum permissible amplicon length: " + str(self.length) + "\n" + \
        "bp position of permissible amplicon start: " + str(self.permissible_start) + "\n" + \
        "Number of guides contained: " + str(self.guidecount) + "\n" + \
        "Required priming subregion: " + str(self.required_start_relative) + "-" + str(self.required_end_relative) + "\n"
        # .start, .end: count of runs along chromosome
        # .fiveprimeabut, .threeprimeabut: count of guides along chromosome that abut good guides
        # .guides: list of good guides in run
        # .length: permissible region distance
        # .guidecount: number of guides in list (=len(self.guides))
        # .genomename: genomename
        # .chromosome
        # .left_outside: the id of the guide to the left of the outside
        # .left_inside: ...
        # .permissible_start: absolute numbering vs chromosome
        # .required_start_absolute: absolute numbering vs chromosome
        # .required_start_relative
        # .permissible_region


