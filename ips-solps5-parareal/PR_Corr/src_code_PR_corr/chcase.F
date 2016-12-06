      subroutine chcase (ind, chtxt)
      use b2mod_types
      integer ind
      character chtxt*(*)
*     ------------------------------------------------------------------
*     CHCASE converts a character variable to either all uppercase or
*     all lowercase.
*!!   The text of this routine is case-sensitive.
*     ------------------------------------------------------------------
      integer k
      intrinsic len, char, ichar
*     ------------------------------------------------------------------
*   ..test installation
*     (protection against errors of case conversion.)
      if (ichar('a').eq.ichar('A')) then
       stop 'chcase--installation error'
      endif
*   ..branch on ind
      if (ind.eq.0) then
*    ..convert chtxt to lowercase
       do k = 1, len(chtxt)
        ichr = ichar(chtxt(k:k))
        if (ichar('A').le.ichr.and.ichr.le.ichar('Z')) then
         chtxt(k:k) = char(ichr+ichar('a')-ichar('A'))
        endif
       enddo
      else if (ind.eq.1) then
*    ..convert chtxt to uppercase
       do k = 1, len(chtxt)
        ichr = ichar(chtxt(k:k))
        if (ichar('a').le.ichr.and.ichr.le.ichar('z')) then
         chtxt(k:k) = char(ichr+ichar('A')-ichar('a'))
        endif
       enddo
      else
       stop 'chcase--faulty argument ind'
      endif
      return
*     ------------------------------------------------------------------
      end
