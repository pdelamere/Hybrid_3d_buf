    MODULE biomath_strings_mod
! Contains procedures (methods) related to string operations
! .. Use Statements ..
      USE biomath_constants_mod
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Public Statements ..
      PUBLIC :: lower_case_char, lower_case_string, qlex, &
        upper_case_char, upper_case_string
! ..
! .. Parameters ..
      CHARACTER (26), PARAMETER :: locase = 'abcdefghijklmnopqrstuvwxyz'
      CHARACTER (26), PARAMETER :: upcase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
! ..
    CONTAINS

!*********************************************************************

      FUNCTION lower_case_char(chr)
!!!-------------------------------------------------------------------
!!!
!!! NAME: translate one character to LOWER case
!!!
!!! DESCRIPTION:
!!!
!!! If CHR is an UPPER case letter 'A..Z', returns the corresponding
!!! lower case letter .
!!! If CHR is not an upper case letter, then returns CHR.
!!!
!!! ARGUMENT:
!!!
!!! CHR --> Character to be translated to LOWER case.
!!!                              CHARACTER*1 CHR
!!!
!!!-------------------------------------------------------------------
! .. Function Return Value ..
        CHARACTER (1) :: lower_case_char
! ..
! .. Scalar Arguments ..
        CHARACTER (1), INTENT (IN) :: chr
! ..
! .. Local Scalars ..
        INTEGER :: ix
        LOGICAL :: qnotuc
! ..
! .. Intrinsic Functions ..
        INTRINSIC INDEX, LGT, LLT
! ..
        qnotuc = (LLT(chr,'A') .OR. LGT(chr,'Z'))

        IF (qnotuc) THEN
! chr is not an UPPER case letter

          lower_case_char = chr

          RETURN

        END IF

        ix = INDEX(upcase,chr)

        IF (ix<=0) THEN
          lower_case_char = chr
        ELSE
! Convert UPPER case to LOWER case
          lower_case_char = locase(ix:ix)
        END IF

        RETURN

      END FUNCTION lower_case_char

!*********************************************************************

      FUNCTION upper_case_char(chr)
!!!-------------------------------------------------------------------
!!!
!!! NAME: translate one character to UPPER case
!!!
!!! DESCRIPTION:
!!!
!!! If CHR is a lower case letter 'a..z' returns the  corresponding
!!! UPPER case letter.
!!! If CHR is not a lower case letter, then returns CHR.
!!!
!!! ARGUMENT:
!!!
!!! CHR --> Character to be translated to UPPER case.
!!!                              CHARACTER*1 CHR
!!!
!!!-------------------------------------------------------------------
! .. Function Return Value ..
        CHARACTER (1) :: upper_case_char
! ..
! .. Scalar Arguments ..
        CHARACTER (1), INTENT (IN) :: chr
! ..
! .. Local Scalars ..
        INTEGER :: ix
        LOGICAL :: qnotuc
! ..
! .. Intrinsic Functions ..
        INTRINSIC INDEX, LGT, LLT
! ..
        qnotuc = (LLT(chr,'a') .OR. LGT(chr,'z'))

        IF (qnotuc) THEN
! chr is not a LOWER case letter

          upper_case_char = chr

          RETURN

        END IF

        ix = INDEX(locase,chr)

        IF (ix<=0) THEN
          upper_case_char = chr
        ELSE
! Convert LOWER case to UPPER case
          upper_case_char = upcase(ix:ix)
        END IF

        RETURN

      END FUNCTION upper_case_char

!*********************************************************************

      SUBROUTINE lower_case_string(str)
!!!-------------------------------------------------------------------
!!!
!!! NAME: LOWer CASE a whole STRING
!!!
!!! DESCRIPTION:
!!!
!!! For each CHR in STR, if CHR is an upper case letter 'A..Z' then it
!!! is replaced by the corresponding lower case letter .
!!! If CHR is not an upper case letter, then nothing happens to CHR.
!!!
!!! ARGUMENT:
!!!
!!! STR <-> String to be translated to lower case.
!!!                             CHARACTER*(*) STR
!!!
!!!-------------------------------------------------------------------
! .. Scalar Arguments ..
        CHARACTER (*), INTENT (INOUT) :: str
! ..
! .. Local Scalars ..
        INTEGER :: ix, j
        LOGICAL :: qnotuc
        CHARACTER (1) :: chr
! ..
! .. Intrinsic Functions ..
        INTRINSIC INDEX, LEN_TRIM, LGT, LLT
! ..
        DO j = 1, LEN_TRIM(str)

          chr = str(j:j)

          qnotuc = (LLT(chr,'A') .OR. LGT(chr,'Z'))

          IF (qnotuc) CYCLE

          ix = INDEX(upcase,chr)

          IF (ix<=0) CYCLE

          str(j:j) = locase(ix:ix)

        END DO

        RETURN

      END SUBROUTINE lower_case_string

!*********************************************************************

      SUBROUTINE upper_case_string(str)
!!!-------------------------------------------------------------------
!!!
!!! NAME: UPper CASE a whole STRING
!!!
!!! DESCRIPTION:
!!!
!!! For each CHR in STR, if CHR is an LOWER case letter 'a..z' then it
!!! is replaced by the corresponding UPPER case letter .
!!! If CHR is not an lower case letter, then nothing happens to CHR.
!!!
!!! ARGUMENT:
!!!
!!! STR <-> String to be translated to upper case.
!!!                             CHARACTER*(*) STR
!!!
!!!-------------------------------------------------------------------
! .. Scalar Arguments ..
        CHARACTER (*), INTENT (INOUT) :: str
! ..
! .. Local Scalars ..
        INTEGER :: ix, j
        LOGICAL :: qnotuc
        CHARACTER (1) :: chr
! ..
! .. Intrinsic Functions ..
        INTRINSIC INDEX, LEN_TRIM, LGT, LLT
! ..
        DO j = 1, LEN_TRIM(str)

          chr = str(j:j)

          qnotuc = (LLT(chr,'a') .OR. LGT(chr,'z'))

          IF (qnotuc) CYCLE

          ix = INDEX(locase,chr)

          IF (ix<=0) CYCLE

          str(j:j) = upcase(ix:ix)

        END DO

        RETURN

      END SUBROUTINE upper_case_string

!*********************************************************************

      FUNCTION qlex(string,qstart,curtyp,cval,lcval,ival,dval,indov)
!------------------------------------------------------------------------
!     Name: Lexical Analyser
!     Description:
!         QLEX is a lexical analyser intended for command languages.
!     QLEX recognizes numbers, names, quoted character strings,
!     delimiters, operators, and two categories of probable error.
!     Ignoring probable errors, the blank character is the universal
!     separator of tokens, and blanks are never part of tokens except
!     for quoted strings.  In particular, a blank token is returned
!     only when it occurs within quotes.  Number and name tokens may be
!     terminated by being followed by a delimiter or operator (or a
!     blank).
!         Complete details on the types of tokens and the process of
!     decomposition are given below in the section titled Details.
!     Use:
!         Each call to QLEX causes it to begin searching for a token
!     starting at the position after the end of the last token returned
!     (if any).  If it finds another token, it will return the value
!     .TRUE., otherwise it will return the value .FALSE.  If QLEX
!     returns .FALSE. then all output variables (except QSTART) are
!     undefined.
!         QSTART is set to .TRUE. on a call to QLEX to indicate that
!     it should begin at the beginning of the string in locating the
!     next token.  Successive calls for more tokens in the same string
!     should ensure that QSTART is .FALSE. on entry to QLEX.  To make
!     this easy, QLEX always changes QSTART to .FALSE. before returning.
!         Since QLEX remembers its position and other things between
!     calls, it should be used only on one string at a time.  Calls
!     which alternate between different strings (with QSTART .FALSE.)
!     will produce unpredictable results.
!     Return Value:
!     QLEX <- .TRUE. if a token has been identified and returned,
!             .FALSE. if there are no tokens remaining in STRING.
!             LOGICAL QLEX
!     Arguments:
!     STRING -> The character string to be decomposed into lexical
!               tokens.
!               CHARACTER*(*) STRING
!     QSTART <-> On input, QSTART should be set to .TRUE. to indicate
!                that the first token in STRING is to be found.  A
!                value of .FALSE. indicates that the token following
!                the one found in the previous call should be returned.
!                NOTE: An erroneous value of .FALSE. can lead to
!                unpredictable results.
!                LOGICAL QSTART
!     CURTYP <- Defined only if QLEX is .TRUE.  The types are
!                 IN -- INteger
!                 RL -- ReaL
!                 ID -- alphanumeric name or IDentifier
!                 QS -- Quoted character String
!                 DL -- DeLimiter
!                 OP -- string of OPerator symbols
!                 OS -- Other String
!                       A string of characters recognized by QLEX but
!                       not fitting into the above categories.
!                       (This is usually an error.)
!                 UC -- Unrecognized Characters
!                       A string of characters not recognized by QLEX,
!                       i.e., characters that do not occur in some
!                       definition above.
!                       (This is usually an error.)
!             CHARACTER*2 CURTYP
!     CVAL <- Defined only if QLEX is .TRUE.  The character value of
!             the token found.  The token will be left adjusted.
!             For safety, CVAL should be as long as STRING.
!             CHARACTER*(*) CVAL
!     LCVAL <- Defined only if QLEX is .TRUE.  The length of the
!              token returned in CVAL.
!              INTEGER LCVAL
!     IVAL <- Defined only if QLEX is .TRUE. and CURTYP.EQ.'IN'
!             (Integer) or CURTYP.EQ.'RL' (Real).
!             The integer value of the token.
!             INTEGER IVAL
!     DVAL <- Defined only if QLEX is .TRUE. and CURTYP.EQ.'IN' or
!             CURTYP.EQ.'RL' (Integer or Real).  The real value of
!             the token stored in a double precision number.
!             DOUBLE PRECISION DVAL
!     INDOV <- Overflow indicator. If CURTYP.EQ.'IN' or CURTYP.EQ.'RL'
!              indicated if the value input is too large to store.
!              The indicators are:
!                  0 -- No overflow (IVAL and DVAL defined)
!                  1 -- INTEGER overflow (DVAL defined)
!                  2 -- REAL overflow (No values defined)
!              INTEGER INDOV
!     Details:
!         Syntax of the tokens is given in extended Backus-Naur form.
!     Square braces, [...], indicate that the contents can appear zero
!     or one times.  The braces, (*...*), indicate that the contents
!     can occur zero or more times.  The construct <BLANK> indicates the
!     blank character, ' '; <CHARACTER> is any legal character (this
!     varies from machine to machine); <SLASH> is the '/' character
!     which is named because it is used as the OR symbol in the
!     description language.
!     -- Integer or Real Strings --
!         Both types can begin with a sign or a digit.  The sign can
!     be separated from the following characters by blanks.  The fol-
!     lowing characters can be only digits for an integer.  A real
!     can begin with a decimal point, '.', followed by a fraction.
!     A real can have an optional exponential part.  This has the form
!     E (or D) followed by a number.
!         Both types of tokens are terminated by being followed by a
!     blank character, an operator character, a delimiter, or a
!     character not defined to QLEX (the latter is usually an error).
!         NOTE: Spaces between the sign and the first character of
!     the remainder of the number are removed in CVAL.
!         <SIGN> ::= +/-
!         <DIGIT> ::= 0/.../9
!         <DIGITST> ::= <DIGIT> (*<DIGIT>*)
!         <INTEGER> ::= [<SIGN>] <DIGITST>
!         <REAL PART> ::= [<SIGN>] [<DIGITST>].<DIGITST> /
!                         [<SIGN>] <DIGITST>.[<DIGITST>]
!         <EXP IND> ::= d/D/e/E
!         <EXP PART> ::= <EXP IND> [<SIGN>] <DIGITST>
!         <REAL> ::= <REAL PART> [<EXP PART>]
!     -- Alphanumeric Name or Identifier --
!         The token must begin with a letter.  Succeeding characters
!     can be letters, digits, or one of the characters
!     '$', '.', or '_' (underscore character).
!         The token is terminated by being followed by a
!     blank character, an operator character, a delimiter, or a
!     character not defined to QLEX (the latter is usually an error).
!         <LETTER> ::= a/.../z/A/.../Z
!         <OTHER ALPHAS> ::= $/./_
!         <EXTENDED ALPHA> ::= <LETTER>/<DIGIT>/<OTHER ALPHAS>
!         <IDENTIFIER> ::= <LETTER>(*<EXTENDED ALPHA>*)
!     -- Quoted String --
!         A quoted string is surrounded by the same quote delimiters,
!     which can be either a single or double quote (' or ").  The
!     appearance of the quote delimiter within a string is indicated
!     by its doubling.  Thus, the string 'don''t' is recognized as
!     the word "don't".  Any characters whatsoever may be in a quoted
!     string, even characters not recognized by QLEX.
!         The quoted string returned by QLEX in CVAL does not include
!     the surrounding single or double quotes. Doubled occurrences of
!     the quoting character within the string are replaced by a single
!     occurrence within CVAL.
!         <QUOTED STRING> ::= '<CHARACTER>(*<CHARACTER>*)' /
!                             "<CHARACTER>(*<CHARACTER>*)"
!     -- Delimiter --
!         A delimiter is a single character from the set specified.
!         <DELIMITER> ::= , / ( / ) / { / } / : / ;
!         In ASCII, square brackets '[' and ']' are added to the set
!     of delimiter characters. These characters do not exist in EBCDIC.
!     -- String of Operators --
!         A string of operators is any contiguous string of operator
!     characters and is terminated by any character that is not an
!     operator.
!         <OPERATOR CHARACTER> ::= + / - / * / <SLASH> / < / > / =
!         <STRING OF OPERATORS>::= <OPERATOR CHARACTER>
!                                       (*<OPERATOR CHARACTER>*)
!     -- String of Recognized Characters --
!         This is a string of characters that are included in some
!     definition above but such that the string does not fall into
!     any of the above defined categories.  The string is terminated by
!     being followed by a blank or an unrecognized character.
!         This category is usually a mistake caused by the ommission of
!     a blank, ending quote, or some other error.
!     -- String of Unrecognized Characters --
!         This is a string of characters that are not included in some
!     definition above.  The string is terminated by being followed by
!     a recognized character.
!     Notes:
!         A plus or minus character (+ or -) is a unary
!     sign if the current call to QLEX has QSTART = .TRUE., or if the
!     preceding token was an operator string, or if the preceding token
!     was one of the delimiters, or if the next character is a digit.
!     In all other cases, the character is considered to be
!     a binary operator.
!----------------------------------------------------------------------
!****
!     BLANK   --  Constant used for comparison.
!****
!****
!     The twelve non-digit classes of characters. Each digit 0-9
!     has itself for a value.
!     CPER    --  The class consisting of the period, '.'
!     CPLU    --  The class consisting of the plus, '+'
!     CMIN    --  The class consisting of the minus, '-'
!     CEXP    --  The class of exponent indicators,
!                 'E', 'e', 'D', and 'd'
!     CALP    --  The class of alphabetic characters not contained
!                 in CEXP
!     COAL    --  The class of other alphabetic characters not
!                 contained in CPLU or CMIN
!     CBLA    --  The class consisting of the space, ' '
!     CDLM    --  The class consisting of the delimiter characters.
!     CQUO    --  The class consisting of the single and double
!                 quotes
!     COPR    --  The class consisting of the operator characters
!     COTH    --  The class consisting of recognized characters
!                 not fitting in one of the above classes
!     CUNR    --  The class of unrecognizable characters strings
!****
!****
!     The set of states of the finite state machine.
!     SNLL    --  This is not a state of the finite state machine.
!                 This is a dummy entry, used when no state transition
!                 will occur.
!     STRT    --  start state; terminated when non-blank encountered.
!     SSGN    --  sign state; a + or - has been encounterd.
!     SDIG    --  digit state; a digit has been encounterd.
!     SPER    --  period state; a . (decimal point) has been encounterd.
!     SFRC    --  fraction state; a digit after the decimal point has
!                 been encounterd.
!     SEEX    --  exponent indicator state; one of the characters
!                 'E', 'e', 'D', or 'd' has been encountered.
!     SSEX    --  exponent sign state; a + or - has been encounterd
!                 after an exponent indicator.
!     SDEX    --  exponent digit state; a digit of the exponent
!                 has been encountered.
!     SID     --  identifier state; an alphabetic character, digit,
!                 or other alphabetic character has been encountered.
!                 (Digits and other alphabetic characters are valid only
!                 if they are not the first character of the token.)
!     SDLM    --  delimiter state; a delimiter character has been
!                 encountered.
!     SSTR    --  string state; a single or double quote has been
!                 encountered.
!     SOPR    --  operator state; an operator character has been
!                 encountered.
!     SOTH    --  other string state; a recognized character not
!                 fitting one of the above states has been encountered.
!     SUNR    --  unrecognizalbe character state; an unrecognizable
!                 character has been encountered.
!****
!****
!     The set of actions to be performed, depending on the state
!     of the finite state machine and the class of the current
!     character.
!     AAIN    --  accumulate integer part of number and string;
!                 set CURTYP = integer
!     AAFR    --  accumulate fractional part of number and string;
!                 set CURTYP = real
!     AAEX    --  accumulate exponent and string;
!     AAST    --  accumulate string;
!     AASR    --  accumulate string;
!                 set CURTYP = real
!     AASO    --  accumulate string;
!                 set CURTYP = other string
!     AAIP    --  accumulate string;
!                 set CURTYP = integer; set SIGN = plus
!     AAIM    --  accumulate string;
!                 set CURTYP = integer; set SIGN = minus
!     AAEP    --  accumulate string;
!                 set ESIGN = positive
!     AAEM    --  accumulate string;
!                 set ESIGN = MINUS
!     AAID    --  accumulate string;
!                 set CURTYP = identifier
!     AANL    --  accumulate nothing;
!                 take no action
!     AASP    --  accumulate string;
!                 set CURTYP = operator string
!     AASU    --  accumulate string;
!                 set CURTYP = unrecognizable character string
!     ATQS    --  accumulate nothing;
!                 set CURTYP = quoted string
!     ACQS    --  check for end of quoted string;
!                 if end of quoted string, return quoted string
!                 else accumulate string
!     ARIN    --  return integer;
!     ARRL    --  return real number;
!     ARID    --  return identifier;
!     ARDL    --  return delimiter;
!     AROP    --  return operator string;
!     AROS    --  return other string;
!     ARUC    --  return unrecognizable character string;
!     ARSO    --  if previous character is blank,
!                 change sign to operator string and return;
!                 else accumulate string and set CURTYP = other string;
!     ARPO    --  if previous character is blank,
!                 change period to other character string;
!****
!****
!     The above are the various token types that may be returned.
!****
!****
!     CURCHR  --  The character currently being processed.
!     PRVCHR  --  The previous character processed.
!     NXTCHR  --  The next character to be processed.
!     STRDLM  --  The delimiter of the current string, if processing
!                 a quoted string.
!     PRVTYP  --  The type of the token returned in any previous call
!                 to this routine.
!****
!****
!     POSN    --  The position in the input string of the character
!                 currently being processed.
!     INTPT   --  The integer part of the number currently being
!                 processed.
!     EXPPT   --  The exponent of the number currently being
!                 processed.
!     SIGN    --  The sign of the integer or real part of the
!                 number currently being processed.
!     ESIGN   --  The sign of the exponent of the number currently
!                 being processed.
!     LENSTR  --  The length of the string being parsed.
!     CHRTYP  --  The type of the current character.
!     NXCTYP  --  The type of the next character.
!     CLASS   --  The class of the current character.
!     STATE   --  The current state of the finite state machine.
!     ACTION  --  Action to be performed.
!     POWER   --  Power of fraction part of current number.
!     EPOW    --  Exponent (absalute value) of fraction part of
!                 current number.
!****
!****
!     NVALUE  --  Dummy argument of statement function ACCUM.
!     ZDIGIT  --  Dummy argument of statement function ACCUM.
!     ACCUM   --  Statement function.
!****
!****
!     IMAX    --  MACHINE CONSTANT, maximum integer value.
!     DMAX    --  MACHINE CONSTANT, maximum double precision value.
!     RMAX    --  MACHINE CONSTANT, maximum real value.
!****
!****
!     Other variables.
!****
!****
!     POWVEC  --  A vector of powers of 10. (1e-1, 1e-2,..., 1e-20)
!                 This is used for an accurate and fast value.
!****
! .. Function Return Value ..
        LOGICAL :: qlex
! ..
! .. Parameters ..
        REAL (dpkind), PARAMETER :: dmax = HUGE(1.0E0_dpkind)
        REAL (dpkind), PARAMETER :: dmax10 = dmax/10.0E0_dpkind
        REAL (dpkind), PARAMETER :: rmax = HUGE(1.0_dpkind)
        REAL (dpkind), PARAMETER :: drmax = rmax
        REAL (dpkind), PARAMETER :: powvec(20) = (/ 1.0E-1_dpkind, &
          1.0E-2_dpkind, 1.0E-3_dpkind, 1.0E-4_dpkind, 1.0E-5_dpkind, &
          1.0E-6_dpkind, 1.0E-7_dpkind, 1.0E-8_dpkind, 1.0E-9_dpkind, &
          1.0E-10_dpkind, 1.0E-11_dpkind, 1.0E-12_dpkind, 1.0E-13_dpkind, &
          1.0E-14_dpkind, 1.0E-15_dpkind, 1.0E-16_dpkind, 1.0E-17_dpkind, &
          1.0E-18_dpkind, 1.0E-19_dpkind, 1.0E-20_dpkind/)
        INTEGER, PARAMETER :: aaem = 10, aaep = 9, aaex = 3, aafr = 2, &
          aaid = 11, aaim = 8, aain = 1, aaip = 7, aanl = 12, aaso = 6, &
          aasp = 13, aasr = 5, aast = 4, aasu = 14, acqs = 16, ardl = 20, &
          arid = 19, arin = 17, arop = 21, aros = 22, arpo = 25, &
          arrl = 18, arso = 24, aruc = 23, atqs = 15, calp = 14, &
          cbla = 16, cdlm = 17, cexp = 13, cmin = 12, coal = 15, &
          copr = 19, coth = 20, cper = 10, cplu = 11, cquo = 18, &
          cunr = 21
        INTEGER, PARAMETER :: imax = HUGE(1)
        INTEGER, PARAMETER :: sdex = 8, sdig = 3, sdlm = 10, seex = 6, &
          sfrc = 5, sid = 9, snll = 0, sopr = 12, soth = 13, sper = 4, &
          ssex = 7, ssgn = 2, sstr = 11, strt = 1, sunr = 14
        CHARACTER (1), PARAMETER :: blank = ' '
        CHARACTER (2), PARAMETER :: typdl = 'DL', typid = 'ID', &
          typin = 'IN', typop = 'OP', typos = 'OS', typqs = 'QS', &
          typrl = 'RL', typuc = 'UC'
        REAL (dpkind), PARAMETER :: dimax = imax
! ..
! .. Scalar Arguments ..
        REAL (dpkind), INTENT (OUT) :: dval
        INTEGER, INTENT (OUT) :: indov, ival, lcval
        LOGICAL, INTENT (INOUT) :: qstart
        CHARACTER (2), INTENT (OUT) :: curtyp
        CHARACTER (*), INTENT (OUT) :: cval
        CHARACTER (*), INTENT (IN) :: string
! ..
! .. Local Scalars ..
        REAL (dpkind) :: exppt, fracpt, intpt, power
        INTEGER :: action, chrtyp, class, epow, esign, i, j, nxctyp, &
          sign, state
        INTEGER, SAVE :: lenstr, posn
        LOGICAL :: qcopy
        CHARACTER (1) :: curchr, nxtchr, strdlm
        CHARACTER (1), SAVE :: prvchr
        CHARACTER (2), SAVE :: prvtyp
! ..
! .. Intrinsic Functions ..
        INTRINSIC HUGE, LEN_TRIM, REAL
! ..
!***
!     If QSTART is TRUE then initialize one time information.
!****
        IF (qstart) THEN
          posn = 1
          prvchr = ' '
          prvtyp = '  '
          lenstr = LEN_TRIM(string) + 1
        END IF
!****
!     Initalize variables for current TOKEN.
!****
        cval = ' '
        lcval = 0

        ival = 0
        dval = zero

        qstart = .FALSE.
        qlex = .FALSE.

! The current character
        curchr = string(posn:posn)

! The character type
        chrtyp = ncqlex(curchr)

        indov = 0

        intpt = zero
        fracpt = zero

        exppt = zero
        epow = 0

        sign = 1
        esign = 1

        state = 1

!-----------------------------------------------------------------------
!     Read one character at a tine, until end of token.
!-----------------------------------------------------------------------

        DO
          IF (posn>lenstr) EXIT

! DMS    If the caracter is not going to be saved to the buffer cval,
!        set qcopy to .FALSE. (states AANL and ATQS)

          qcopy = .TRUE.

          IF (chrtyp<=9) THEN
! A digit (0..9) was read)
            class = 1
          ELSE
! Anything else other than digit
            class = chrtyp - 8
          END IF

          action = get_action(class,state)

          state = next_state(class,state)

! Get NEXT character
          IF (posn<lenstr) THEN
            nxtchr = string(posn+1:posn+1)
          ELSE
            nxtchr = ' '
          END IF

! Get next caharacter's type
          nxctyp = ncqlex(nxtchr)

          IF ((chrtyp==cplu) .OR. (chrtyp==cmin)) THEN
            IF (prvchr==' ') THEN
              IF (class==(cplu-8)) THEN
                action = aaip
              ELSE
                action = aaim
              END IF
            END IF
          END IF

!-----------------------------------------------------------------------
!     Perform action specified in action table.
!-----------------------------------------------------------------------

          SELECT CASE (action)

          CASE (aain)
!****
!     AAIN    --  accumulate integer part of number and string;
!                 set CURTYP = integer
!****
            intpt = accum(intpt,chrtyp)

            curtyp = typin

          CASE (aafr)
!****
!     AAFR    --  accumulate fractional part of number and string;
!                 set CURTYP = real
!****
            epow = epow + 1

            IF (epow<=20) THEN
              power = powvec(epow)
            ELSE
              power = power*0.1E0_dpkind
            END IF

            fracpt = fracpt + (REAL(chrtyp,kind=dpkind)*power)

            curtyp = typrl

          CASE (aaex)
!****
!     AAEX    --  accumulate exponent and string;
!****
            exppt = accum(exppt,chrtyp)

          CASE (aast)
!****
!     AAST    --  accumulate string;
!****

          CASE (aasr)
!****
!     AASR    --  accumulate string;
!                 set CURTYP = real
!****
            curtyp = typrl

          CASE (aaso)
!****
!     AASO    --  accumulate string;
!                 set CURTYP = other string
!****
            curtyp = typos

          CASE (aaip)
!****
!     AAIP    --  accumulate string;
!                 set CURTYP = integer; set SIGN = plus
!****
            curtyp = typin
            sign = 1

          CASE (aaim)
!****
!     AAIM    --  accumulate string;
!                 set CURTYP = integer; set SIGN = minus
!****
            curtyp = typin
            sign = -1

          CASE (aaep)
!****
!     AAEP    --  accumulate string;
!                 set ESIGN = positive
!****
            esign = 1

          CASE (aaem)
!****
!     AAEM    --  accumulate string;
!                 set ESIGN = MINUS
!****
            esign = -1

          CASE (aaid)
!****
!     AAID    --  accumulate string;
!                 set CURTYP = identifier
!****
            curtyp = typid

          CASE (aanl)
!****
!     AANL    --  accumulate nothing;
!                 take no action
!****
            qcopy = .FALSE.

          CASE (aasp)
!****
!     AASP    --  accumulate string;
!                 set CURTYP = operator string
!****
            curtyp = typop

          CASE (aasu)
!****
!     AASU    --  accumulate string;
!                 set CURTYP = unrecognizable character string
!****
            curtyp = typuc

          CASE (atqs)
!****
!     ATQS    --  accumulate nothing;
!                 set CURTYP = quoted string
!****
            strdlm = curchr
            curtyp = typqs
            qcopy = .FALSE.

          CASE (acqs)
!****
!     ACQS    --  check for end of quoted string;
!                 if end of quoted string, return quoted string
!                 else accumulate string
!****
            IF (curchr==strdlm) THEN
              IF (string(posn+1:posn+1)==strdlm) THEN
                posn = posn + 1
              ELSE
                posn = posn + 1
                prvtyp = curtyp
                qlex = .TRUE.

                RETURN
              END IF
            END IF

          CASE (arin)
!****
!     ARIN    --  return integer;
!****
            CALL set_numeric_values

            RETURN

          CASE (arrl)
!****
!     ARRL    --  return real number;
!****
            CALL set_numeric_values

            RETURN

          CASE (arid)
!****
!     ARID    --  return identifier;
!****
            prvtyp = curtyp
            qlex = .TRUE.

            RETURN

          CASE (ardl)
!****
!     ARDL    --  return delimiter;
!****
            lcval = lcval + 1
            cval(lcval:lcval) = curchr
            curtyp = typdl
            posn = posn + 1
            prvtyp = curtyp
            qlex = .TRUE.

            RETURN

          CASE (arop)
!****
!     AROP    --  return operator string;
!****
            prvtyp = curtyp
            qlex = .TRUE.

            RETURN

          CASE (aros)
!****
!     AROS    --  return other string;
!****
            prvtyp = curtyp
            qlex = .TRUE.

            RETURN

          CASE (aruc)
!****
!     ARUC    --  return unrecognizable character string;
!****
            prvtyp = curtyp
            qlex = .TRUE.

            RETURN

          CASE (arso)
!****
!     ARSO    --  if previous character is blank,
!                 change sign to operator string and return;
!                 else accumulate string and set CURTYP = other string;
!****
            IF (prvchr==blank) THEN
              posn = posn - 1
              curtyp = typop
              prvtyp = curtyp
              qlex = .TRUE.
              RETURN

            ELSE
              curtyp = typos
            END IF

          CASE (arpo)
!****
!     ARPO    --  if previous character is blank,
!                 change period to other character string and return;
!                 else return real;
!****
            IF (cval(1:1)=='.') THEN
              curtyp = typos
              prvtyp = curtyp
              qlex = .TRUE.

              RETURN

            ELSE
              CALL set_numeric_values

              RETURN

            END IF

          END SELECT

          IF (qcopy) THEN
            lcval = lcval + 1
            cval(lcval:lcval) = curchr
          END IF

!        If this is the first character encountered and it is a plus,
!        '+' or minus '-' and the previous character was a separator
!        or the next character is a digit then set the SIGN state.

          IF ((chrtyp==cplu) .OR. (chrtyp==cmin)) THEN

            IF (prvchr==' ') THEN
              state = ssgn
            END IF
          END IF

          prvchr = curchr
          curchr = nxtchr
          chrtyp = nxctyp

          posn = posn + 1

        END DO

        IF (curtyp==typqs) THEN
          curtyp = typos
          prvtyp = curtyp
          qlex = .TRUE.
        END IF

        IF (state==ssgn) THEN
          curtyp = typop
          prvtyp = curtyp
          qlex = .TRUE.
        END IF
!     ..
!     .. Return Statement ..
        RETURN

      CONTAINS

!.....................................................................

        FUNCTION accum(nvalue,zdigit)
!     This function is used to accumulate the integer part
!     or the exponent of the number currently being processed.
! ..
! .. Function Return Value ..
          REAL (dpkind) :: accum
! ..
! .. Parameters ..
          REAL (dpkind), PARAMETER :: ten = 10.0_dpkind
! ..
! .. Scalar Arguments ..
          REAL (dpkind), INTENT (IN) :: nvalue
          INTEGER, INTENT (IN) :: zdigit
! ..
          accum = (ten*nvalue) + zdigit

        END FUNCTION accum

!.....................................................................

        FUNCTION get_action(class,state)
! STRT - start state
! SSGN - sign state
! SDIG - digit state
! SPER - period state
! SFRC - fraction state
! SEEX - exponent indicator state
! SSEX - exponent sign state
! SDEX - exponent digit state
! SID  - identifier state
! SDLM - delimiter state
! SSTR - string state
! SOPR - operator state
! SOTH - other character state
! SUNR - unrecognizable character state
! ..
! .. Function Return Value ..
          INTEGER :: get_action
! ..
! .. Parameters ..
          INTEGER, PARAMETER :: actntb(13,14) = RESHAPE((/aain,aasr,aasp, &
            aasp,aaid,aaid,aaso,aanl,ardl,atqs,aasp,aaso,aasu,aain,aasr, &
            arso,arso,arso,arso,arso,aanl,arso,arso,arso,arso,arso,aain, &
            aasr,arin,arin,aasr,aaso,aaso,arin,arin,aaso,arin,aaso,arin, &
            aafr,aaso,arrl,arrl,aast,aaso,aaso,arpo,arrl,aaso,arrl,aaso, &
            arrl,aafr,aaso,arrl,arrl,aasr,aaso,aaso,arrl,arrl,aaso,arrl, &
            aaso,arrl,aaex,aaso,aaep,aaem,aaso,aaso,aaso,aros,aaso,aaso, &
            aaso,aaso,aros,aaex,aaso,aaso,aaso,aaso,aaso,aaso,aros,aaso, &
            aaso,aaso,aaso,aros,aaex,aaso,arrl,arrl,aaso,aaso,aaso,arrl, &
            arrl,aaso,arrl,aaso,arrl,aast,aast,arid,arid,aast,aast,aast, &
            arid,arid,aaso,arid,aaso,arid,aanl,aanl,aanl,aanl,aanl,aanl, &
            aanl,aanl,aanl,aanl,aanl,aanl,aanl,aast,aast,aast,aast,aast, &
            aast,aast,aast,aast,acqs,aast,aast,aast,arop,arop,aast,aast, &
            arop,arop,arop,arop,arop,arop,aast,arop,arop,aast,aast,aast, &
            aast,aast,aast,aast,aros,aast,aast,aast,aast,aros,aruc,aruc, &
            aruc,aruc,aruc,aruc,aruc,aruc,aruc,aruc,aruc,aruc,aast/), &
            shape=(/13,14/),order=(/1,2/))
! ..
! .. Scalar Arguments ..
          INTEGER, INTENT (IN) :: class, state
! ..
! .. Intrinsic Functions ..
          INTRINSIC RESHAPE
! ..
          get_action = actntb(class,state)

          RETURN

        END FUNCTION get_action

!.....................................................................

        FUNCTION ncqlex(char)
!****
!     This function looks up the class of a given character
!     in the CHRSET array.
!****
! 8 uprintable characters
! 8 uprintable characters
! 8 uprintable characters
! 8 uprintable characters
!            BLANK  !     "     #    $     %     &     '
!            (     )     *     +     ,     -     .     /
!            0     1      2     3     4     5     6     7
!             8     9     :     ;     <     =     >     ?
!             @     A     B     C     D     E     F     G
!             H     I     J     K     L     M     N     O
!             P     Q     R     S     T     U     V     W
!             X     Y     Z     [     \     ]     ^     _
!             `     a     b     c     d     e     f     g
!             h     i     j     k     l     m     n     o
!             p     q     r     s     t     u     v     w
! ..
! .. Function Return Value ..
          INTEGER :: ncqlex
! ..
! .. Parameters ..
          INTEGER, PARAMETER :: chrset(0:127) = (/ cunr, cunr, cunr, &
            cunr, cunr, cunr, cunr, cunr, cunr, cunr, cunr, cunr, cunr, &
            cunr, cunr, cunr, cunr, cunr, cunr, cunr, cunr, cunr, cunr, &
            cunr, cunr, cunr, cunr, cunr, cunr, cunr, cunr, cunr, cbla, &
            cunr, cquo, cunr, coal, cunr, cunr, cquo, cdlm, cdlm, copr, &
            cplu, cdlm, cmin, cper, copr, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, &
            cdlm, cdlm, copr, copr, copr, cunr, cunr, calp, calp, calp, &
            cexp, cexp, calp, calp, calp, calp, calp, calp, calp, calp, &
            calp, calp, calp, calp, calp, calp, calp, calp, calp, calp, &
            calp, calp, calp, cdlm, cunr, cdlm, cunr, coal, cunr, calp, &
            calp, calp, cexp, cexp, calp, calp, calp, calp, calp, calp, &
            calp, calp, calp, calp, calp, calp, calp, calp, calp, calp, &
            calp, calp, calp, calp, calp, cdlm, cunr, cdlm, cunr, cunr/)
! ..
! .. Scalar Arguments ..
          CHARACTER (1), INTENT (IN) :: char
! ..
! .. Intrinsic Functions ..
          INTRINSIC IACHAR
! ..
!             x     y     z     {     |     }     ~     DEL
          ncqlex = chrset(IACHAR(char))

          RETURN
        END FUNCTION ncqlex

!.....................................................................

        FUNCTION next_state(class,state)
!****
!     Transition Table.
!     Row indices are the character classes
!     Column index is the current state
!     Table entry is next state given current state and class of
!     the next character seen
!****
!!! Row indices of trantab
!!!   dig  cper cplu cmin cexp calp coal cbla cdlm cquo copr coth cunr
! STRT - start state
! SSGN - sign state
! SDIG - digit state
! SPER - period state
! SFRC - fraction state
! SEEX - exponent indicator state
! SSEX - exponent sign state
! SDEX - exponent digit state
! SID  - identifier state
! SDLM - delimiter state
! SSTR - string state
! SOPR - operator state
! SOTH - other character state
! SUNR - unrecognizable character state
! ..
! .. Function Return Value ..
          INTEGER :: next_state
! ..
! .. Parameters ..
          INTEGER, PARAMETER :: trantb(13,14) = RESHAPE((/sdig,sper,sopr, &
            sopr,sid,sid,soth,strt,snll,sstr,sopr,soth,sunr,sdig,sper, &
            snll,snll,snll,snll,snll,ssgn,snll,snll,snll,snll,snll,sdig, &
            sper,snll,snll,seex,soth,soth,snll,snll,soth,snll,soth,snll, &
            sfrc,soth,snll,snll,seex,soth,soth,snll,snll,soth,snll,soth, &
            snll,sfrc,soth,snll,snll,seex,soth,soth,snll,snll,soth,snll, &
            soth,snll,sdex,soth,ssex,ssex,soth,soth,soth,snll,soth,soth, &
            soth,soth,snll,sdex,soth,soth,soth,soth,soth,soth,snll,soth, &
            soth,soth,soth,snll,sdex,soth,snll,snll,soth,soth,soth,snll, &
            snll,soth,snll,soth,snll,sid,sid,snll,snll,sid,sid,sid,snll, &
            snll,soth,snll,soth,snll,snll,snll,snll,snll,snll,snll,snll, &
            snll,snll,snll,snll,snll,snll,sstr,sstr,sstr,sstr,sstr,sstr, &
            sstr,sstr,sstr,sstr,sstr,sstr,sstr,snll,snll,sopr,sopr,snll, &
            snll,snll,snll,snll,snll,sopr,snll,snll,soth,soth,soth,soth, &
            soth,soth,soth,snll,soth,soth,soth,soth,snll,snll,snll,snll, &
            snll,snll,snll,snll,snll,snll,snll,snll,snll,sunr/), &
            shape=(/13,14/))
! ..
! .. Scalar Arguments ..
          INTEGER, INTENT (IN) :: class, state
! ..
! .. Intrinsic Functions ..
          INTRINSIC RESHAPE
! ..
!!! Row indices of trantab
!!!   dig  cper cplu cmin cexp calp coal cbla cdlm cquo copr coth cunr
          next_state = trantb(class,state)

        END FUNCTION next_state

!.....................................................................

        SUBROUTINE set_numeric_values
! .. Intrinsic Functions ..
          INTRINSIC INT
! ..
          IF (curtyp==typin) THEN
! INTEGER value

            dval = sign*intpt

            IF (intpt<=dimax) THEN
              indov = 0
              ival = sign*INT(intpt)
            ELSE IF (intpt<=drmax) THEN
              indov = 1
              ival = 0
            ELSE
              indov = 2
              ival = 0
            END IF
          ELSE IF (curtyp==typrl) THEN
! REAL value

            dval = intpt + fracpt

            IF (exppt>=REAL(imax,kind=dpkind)) THEN
              indov = 3
              ival = 0
              dval = zero
              prvtyp = curtyp
              qlex = .TRUE.

              RETURN

            END IF

            i = INT(exppt)
            DO j = 1, i
              IF (esign==1) THEN
                IF (dval<=dmax10) THEN
                  dval = dval*ten
                ELSE
                  indov = 3
                  ival = 0
                  dval = zero
                  prvtyp = curtyp
                  qlex = .TRUE.

                  RETURN

                END IF
              ELSE
                dval = dval*tenth
              END IF
            END DO

            indov = 2

            IF (dval<=REAL(rmax,kind=dpkind)) THEN
              indov = 1

              IF (dval<=dimax) THEN
                indov = 0
                ival = INT(dval)
              END IF
            END IF

! Change the sign if the string starts with a '-'

            IF (cval(1:1)=='-') THEN
              ival = -ival
              dval = -dval
            END IF
          ELSE
            ival = 0
            dval = zero
          END IF

!-----------------------------------------------------------------------
!     Return with a TOKEN.
!-----------------------------------------------------------------------
          prvtyp = curtyp
          qlex = .TRUE.

        END SUBROUTINE set_numeric_values

!.....................................................................

      END FUNCTION qlex

!*********************************************************************

    END MODULE biomath_strings_mod
