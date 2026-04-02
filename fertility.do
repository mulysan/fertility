 import delimited "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility-rate-by-age-group.csv", encoding(UTF-8)clear 
 keep if entity == "Israel"
 
 reshape long f, i(entity code year) j(age)
 
 replace f = f*5/1000
 rename f Births
 gen birth_year = year -age
 keep if mod(birth_year, 5) ==1
 
 *bys birth_year: keep if _N == 9
 *collapse (sum) f, by(birth_year)
 
 keep if age > 20 & age < 45
 levelsof age, local(years)
local plots ""
local labels ""
local i = 1
foreach y of local years {
    local plots "`plots' (line Births birth_year if age==`y')"
	local y1 = `y' - 2
	local y2 = `y' + 2

    local labels "`labels' label(`i' "Age `y1' to `y2'")"
    local i = `i' + 1
}

label var birth_year "Birth year"
twoway `plots', legend(`labels') ///
       title("Births by Birth Year and Age: Israel") ///
       xlabel(, angle(45))
	   
	   
	   graph export "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility_birthyear_age.png", replace

	   * Chile
	   
	    import delimited "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility-rate-by-age-group.csv", encoding(UTF-8)clear 
 keep if entity == "Chile"
 
 reshape long f, i(entity code year) j(age)
 
 replace f = f*5/1000
 rename f Births
 gen birth_year = year -age
 keep if mod(birth_year, 5) ==1
 
 *bys birth_year: keep if _N == 9
 *collapse (sum) f, by(birth_year)
 
 keep if age > 20 & age < 45
 levelsof age, local(years)
local plots ""
local labels ""
local i = 1
foreach y of local years {
    local plots "`plots' (line Births birth_year if age==`y')"
	local y1 = `y' - 2
	local y2 = `y' + 2

    local labels "`labels' label(`i' "Age `y1' to `y2'")"
    local i = `i' + 1
}

label var birth_year "Birth year"
twoway `plots', legend(`labels') ///
       title("Births by Birth Year and Age: Chile") ///
       xlabel(, angle(45))
	   
	   
	   graph export "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility_birthyear_age_chile.png", replace

	   
	   
bys birth_year: egen total_births = sum(Births)
bys birth_year: gen obs_sum = _N
gen share51 = Births/total_births if birth_year == 1951
bys age: egen share511 = min(share51)


gen Births1 = total_births * share511

collapse (sum) Births Births1 (count) obs = Births (mean) obs_sum , by(year)
*keep if year >= 1993 & year <= 2003

keep if obs_sum == 5

twoway (line Births year ) (line Births1 year ), ///
       title("Actual and Counterfactual TFR by year") legend( label(1 "Actual TFR") label(2 "Counterfactual TFR") ) ///
       xlabel(, angle(45))
	   
	   
	   graph export "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility_counterfactual.png", replace

*
	   
	 import delimited "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility-rate-by-age-group.csv", encoding(UTF-8)clear 
 keep if entity == "Chile"
 
 reshape long f, i(entity code year) j(age)
 
 replace f = f*5/1000
 rename f Births
 gen birth_year = year -age
 keep if mod(birth_year, 5) ==1

 keep if age > 20 & age < 40
 
 bys birth_year: egen total_births = sum(Births)
bys birth_year: gen obs_sum = _N
gen share51 = Births/total_births if birth_year == 1951
bys age: egen share511 = min(share51)


gen Births1 = total_births * share511

collapse (sum) Births Births1 (count) obs = Births (mean) obs_sum , by(year)
*keep if year >= 1993 & year <= 2003

keep if obs_sum == 4

twoway (line Births year ) (line Births1 year ), ///
       title("Actual and Counterfactual TFR by year") legend( label(1 "Actual TFR") label(2 "Counterfactual TFR") ) ///
       xlabel(, angle(45))
	   
	   
	   graph export "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility_counterfactual_20_40.png", replace

	   *
	   
	 import delimited "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility-rate-by-age-group.csv", encoding(UTF-8)clear 
 keep if entity == "Israel"
 
 reshape long f, i(entity code year) j(age)
 
 replace f = f*5/1000
 rename f Births
 gen birth_year = year -age
 *keep if mod(birth_year, 5) ==1

 keep if age >= 20 & age <= 39
 
 bys birth_year: egen total_births = sum(Births)
bys birth_year: gen obs_sum = _N
gen share51 = Births/total_births if birth_year == 1951
bys age: egen share511 = min(share51)


gen Births1 = total_births * share511

collapse (sum) Births Births1 (count) obs = Births (mean) obs_sum , by(year)
keep if year >= 1970

twoway (line Births year ) (line Births1 year if obs_sum == 4), ///
       title("Actual and Counterfactual TFR by year") legend( label(1 "Actual TFR") label(2 "Counterfactual TFR") ) ///
       xlabel(, angle(45))
	   
	   
	   graph export "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility_counterfactual_20_40_israel.png", replace

	  *
	   
	 import delimited "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility-rate-by-age-group.csv", encoding(UTF-8)clear 
 keep if entity == "Israel"
 
 reshape long f, i(entity code year) j(age)
 
 replace f = f*5/1000
 rename f Births

 
 collapse (sum) Births   , by(year)

 twoway (line Births year ), ///
       title(" TFR by year")  ///
       xlabel(, angle(45))
	   
	   
	   graph export "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility_israel.png", replace

	   
	     *
	   
	 import delimited "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/fertility-rate-by-age-group.csv", encoding(UTF-8)clear 
 keep if entity == "Israel"
 
 reshape long f, i(entity code year) j(age)
 
 replace f = f*5/1000
 rename f Births
 
 
 bys year: egen total_births = sum(Births)
 gen share = Births/total_births 

gen age_w = age * share
 
 collapse (sum) age_w   , by(year)

 label var age_w "Age"
 
 twoway (line age_w year ), ///
       title("Averge age at birth: Israel")  ///
       xlabel(, angle(45))
	   
	   
	   graph export "/Users/shmuelsan/Dropbox/HUJI/misc/twitter/fertility/average_age_israel.png", replace

	   
	   
