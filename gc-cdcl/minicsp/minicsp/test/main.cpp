/*************************************************************************
minicsp

Copyright 2010--2014 George Katsirelos

Minicsp is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Minicsp is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with minicsp.  If not, see <http://www.gnu.org/licenses/>.

*************************************************************************/

void bt_test();
void dom_test();
void schedule_test();
void pb_test();
void le_test();
void lin_test();
void element_test();
void regular_test();
void alldiff_test();
void nvalue_test();
void set_test();
void lex_test();

int main()
{
  bt_test();
  dom_test();
  schedule_test();
  pb_test();
  le_test();
  lin_test();
  element_test();
  regular_test();
  alldiff_test();
  nvalue_test();
  set_test();
  lex_test();
  return 0;
}
