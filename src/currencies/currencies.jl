### Based off of Ito.jl
type NullCurrency <: AbstractCurrency end

immutable Currency{S <: AbstractString, S2 <: AbstractString, S3 <: AbstractString, I <: Integer} <: AbstractCurrency
	name::S
	code::S
	numeric::I
	symbol::S2
	fractionSymbol::S3
	fractionsPerUnit::I
	rounding::Function
	formatString::S
end

# Data from http://fx.sauder.ubc.ca/currency_table.html
# and http://www.thefinancials.com/vortex/CurrencyFormats.html

list_currencies=[

# Africa
("South-African rand", "ZAR", 710, "R", "", 100, identity, "%3% %1\$.2f"),

# America
("Argentinian peso", "ARS", 32, "", "", 100, identity, "%2% %1\$.2f"),
("Brazilian real", "BRL", 986, "R\$", "", 100, identity, "%3% %1\$.2f"),
("Canadian dollar", "CAD", 124, "Can\$", "", 100, identity, "%3% %1\$.2f"),
("Chilean peso", "CLP", 152, "Ch\$", "", 100, identity, "%3% %1\$.0f"),
("Colombian peso", "COP", 170, "Col\$", "", 100, identity, "%3% %1\$.2f"),
("Mexican peso", "MXN", 484, "Mex\$", "", 100, identity, "%3% %1\$.2f"),
("Peruvian nuevo sol", "PEN", 604, "S/.", "", 100, identity, "%3% %1\$.2f"),
("Peruvian inti", "PEI", 998, "I/.", "", 100, identity, "%3% %1\$.2f"),
("Peruvian sol", "PEH", 999, "S./", "", 100, identity, "%3% %1\$.2f"),
("Trinidad & Tobago dollar", "TTD", 780, "TTD\$", "", 100, identity, "%3% %1\$.2f"),
("U.S. dollar", "USD", 840, "\$", "¢", 100, identity, "%3% %1\$.2f"),
("Venezuelan bolivar", "VEB", 862, "Bs", "", 100, identity, "%3% %1\$.2f"),

# Asia
("Bangladesh taka", "BDT", 50, "Bt", "", 100, identity, "%3% %1\$.2f"),
("Chinese yuan", "CNY", 156, "Y", "", 100, identity, "%3% %1\$.2f"),
("Hong Kong dollar", "HKD", 344, "HK\$", "", 100, identity, "%3% %1\$.2f"),
("Israeli shekel", "ILS", 376, "NIS", "", 100, identity, "%1\$.2f %3%"),
("Indian rupee", "INR", 356, "Rs", "", 100, identity, "%3% %1\$.2f"),
("Iraqi dinar", "IQD", 368, "ID", "", 1000, identity, "%2% %1\$.3f"),
("Iranian rial", "IRR", 364, "Rls", "", 1, identity, "%3% %1\$.2f"),
("Japanese yen", "JPY", 392, "¥", "", 100, identity, "%3% %1\$.0f"),
("South-Korean won", "KRW", 410, "W", "", 100, identity, "%3% %1\$.0f"),
("Kuwaiti dinar", "KWD", 414, "KD", "", 1000, identity, "%3% %1\$.3f"),
("Nepal rupee", "NPR", 524, "NRs", "", 100, identity, "%3% %1\$.2f"),
("Pakistani rupee", "PKR", 586, "Rs", "", 100, identity, "%3% %1\$.2f"),
("Saudi riyal", "SAR", 682, "SRls", "", 100, identity, "%3% %1\$.2f"),
("Singapore dollar", "SGD", 702, "S\$", "", 100, identity, "%3% %1\$.2f"),
("Thai baht", "THB", 764, "Bht", "", 100, identity, "%1\$.2f %3%"),
("Taiwan dollar", "TWD", 901, "NT\$", "", 100, identity, "%3% %1\$.2f"),

# Europe
("Bulgarian lev", "BGL", 100, "lv", "", 100, identity, "%1\$.2f %3%"),
("Belarussian ruble", "BYR", 974, "BR", "", 1, identity, "%2% %1\$.0f"),
("Swiss franc", "CHF", 756, "SwF", "", 100, identity, "%3% %1\$.2f"),
("Czech koruna", "CZK", 203, "Kc", "", 100, identity, "%1\$.2f %3%"),
("Danish krone", "DKK", 208, "Dkr", "", 100, identity, "%3% %1\$.2f"),
("European Euro", "EUR", 978, "", "", 100, x->round(x,2), "%2% %1\$.2f"),
("British pound sterling", "GBP", 826,  "£", "p", 100, identity, "%3% %1\$.2f"),
("Hungarian forint", "HUF", 348, "Ft", "", 1, identity, "%1\$.0f %3%"),
("Iceland krona", "ISK", 352, "IKr", "", 100, identity, "%1\$.2f %3%"),
("Lithuanian litas", "LTL", 440, "Lt", "", 100, identity, "%1\$.2f %3%"),
("Norwegian krone", "NOK", 578, "NKr", "", 100, identity, "%3% %1\$.2f"),
("Polish zloty", "PLN", 985, "zl", "", 100, identity, "%1\$.2f %3%"),
("Romanian leu", "ROL", 642, "L", "", 100, identity, "%1\$.2f %3%"),
("Romanian new leu", "RON", 946, "L", "", 100, identity, "%1\$.2f %3%"),
("Swedish krona", "SEK", 752, "kr", "", 100, identity, "%1\$.2f %3%"),
("Turkish lira", "TRL", 792, "TL", "", 100, identity, "%1\$.0f %3%"),
("New Turkish lira", "TRY", 949, "YTL", "", 100, identity, "%1\$.2f %3%"),

# Europe deprecated
("Austrian shilling", "ATS", 40, "", "", 100, identity, "%2% %1\$.2f"),
("Belgian franc", "BEF", 56, "", "", 1, identity, "%2% %1\$.0f"),
("Cyprus pound", "CYP", 196, "£", "", 100, identity, "%3% %1\$.2f"),
("Deutsche mark", "DEM", 276, "DM", "", 100, identity, "%1\$.2f %3%"),
("Estonian kroon", "EEK", 233, "KR", "", 100, identity, "%1\$.2f %2%"),
("Spanish peseta", "ESP", 724, "Pta", "", 100, identity, "%1\$.0f %3%"),
("Finnish markka", "FIM", 246, "mk", "", 100, identity, "%1\$.2f %3%"),
("French franc", "FRF", 250, "", "", 100, identity, "%1\$.2f %2%"),
("Greek drachma", "GRD", 300, "", "", 100, identity, "%1\$.2f %2%"),
("Irish punt", "IEP", 372, "", "", 100, identity, "%2% %1\$.2f"),
("Italian lira", "ITL", 380, "L", "", 1, identity, "%3% %1\$.0f"),
("Luxembourg franc", "LUF", 442, "F", "", 100, identity, "%1\$.0f %3%"),
("Latvian lat", "LVL", 428, "Ls", "", 100, identity, "%3% %1\$.2f"),
("Maltese lira", "MTL", 470, "Lm", "", 100, identity, "%3% %1\$.2f"),
("Dutch guilder", "NLG", 528, "f", "", 100, identity, "%3% %1\$.2f"),
("Portuguese escudo", "PTE", 620, "Esc", "", 100, identity, "%1\$.0f %3%"),
("Slovenian tolar", "SIT", 705, "SlT", "", 100, identity, "%1\$.2f %3%"),
("Slovak koruna", "SKK", 703, "Sk", "", 100, identity, "%1\$.2f %3%"),

# Oceania
("Australian dollar", "AUD", 36, "A\$", "", 100, identity, "%3% %1\$.2f"),
("New Zealand dollar", "NZD", 554, "NZ\$", "", 100, identity, "%3% %1\$.2f")
]

list_deprecated=Dict(
"ATS"=>"EUR",
"BEF"=>"EUR",
"CYP"=>"EUR",
"DEM"=>"EUR",
"EEK"=>"EUR",
"ESP"=>"EUR",
"FIM"=>"EUR",
"FRF"=>"EUR",
"GRD"=>"EUR",
"IEP"=>"EUR",
"ITL"=>"EUR",
"LUF"=>"EUR",
"LVL"=>"EUR",
"MTL"=>"EUR",
"NLG"=>"EUR",
"PTE"=>"EUR",
"SIT"=>"EUR",
"SKK"=>"EUR",
"TRL"=>"TRY",
"ROL"=>"RON",
"PEH"=>"PEI",
"PEI"=>"PEN"
)

# Codegen function
for currency in list_currencies
	@eval ($(symbol("$(currency[2])"*"Currency")))()=Currency(($currency)...)
end
