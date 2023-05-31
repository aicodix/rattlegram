/*
Some constants

Constants below are the result of truncating the output of "bc -l" computed with "scale=100".

Copyright 2018 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

namespace DSP {

template <typename T>
struct Const
{
	static constexpr T EighthPi()
	{
		// a(1)/2
		return .3926990816987241548078304229099378605246461749218882276218680740384L;
	}
	static constexpr T FourthPi()
	{
		// a(1)
		return .7853981633974483096156608458198757210492923498437764552437361480769L;
	}
	static constexpr T HalfPi()
	{
		// 2*a(1)
		return 1.570796326794896619231321691639751442098584699687552910487472296153L;
	}
	static constexpr T Pi()
	{
		// 4*a(1)
		return 3.141592653589793238462643383279502884197169399375105820974944592307L;
	}
	static constexpr T TwoPi()
	{
		// 8*a(1)
		return 6.283185307179586476925286766559005768394338798750211641949889184615L;
	}
	static constexpr T FourPi()
	{
		// 16*a(1)
		return 12.56637061435917295385057353311801153678867759750042328389977836923L;
	}
	static constexpr T SqrtPi()
	{
		// sqrt(4*a(1))
		return 1.772453850905516027298167483341145182797549456122387128213807789852L;
	}
	static constexpr T SqrtTwoPi()
	{
		// sqrt(8*a(1))
		return 2.506628274631000502415765284811045253006986740609938316629923576342L;
	}
	static constexpr T InvSqrtTwo()
	{
		// sqrt(2)/2
		return .7071067811865475244008443621048490392848359376884740365883398689953L;
	}
};

}

