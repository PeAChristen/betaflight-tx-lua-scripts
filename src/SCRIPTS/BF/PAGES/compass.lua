local template = assert(loadScript(radio.template))()
local margin = template.margin
local indent = template.indent
local lineSpacing = template.lineSpacing
local tableSpacing = template.tableSpacing
local sp = template.listSpacing.field
local yMinLim = radio.yMinLimit
local x = margin
local y = yMinLim - lineSpacing
local inc = { x = function(val) x = x + val return x end, y = function(val) y = y + val return y end }
local labels = {}
local fields = {}

local path_to_wmm = "/SCRIPTS/BF/COMPASS/WMM.COF"
local wmm_data = {}
local gps_data = {lat=0, long 0, alt=0, sats=0}
local position_data = {lat=0, long=0, alt=0}
local declination = 0

local MSP_RAW_GPS = 106  -- out message: Fix, numsat, lat, lon, alt, speed, ground course
local GPS_TELEM_NAME = "GPS"
local GPS_SATS_TELEM_NAME = "Sats"
local GPS_ALT_TELEM_NAME = "Alt"


-- TODO: Make user to choose from gps data or manual entered data or skip to enter magnetic decliantion directly
-- TODO: Make gps poll loop and show num of gps sats of requierd and be able to escape back to above choice

-- Load WMM.COF, if no file - display info where to put WMM.COF file /SCRIPTS/BF/COMPASS
-- If in gps_data mode and not enough satelites, dispaly wait message, use betaflight min sats value
-- If not enough satelites, user should be able to escape back to manualy enter lat, long and alt
-- Calculate magnetic declination when gps or user data is avalible
-- Let user evaluate Lat, Long, Alt and calculated magnetic declination
-- Let user save the new calculated magnetic delcination value to FC or manualy enterd value

labels[#labels + 1] = { t = "Compass Settings",    x = x, y = inc.y(lineSpacing) }
fields[#fields + 1] = { t = "Magnetic Declination",     x = x + indent, y = inc.y(lineSpacing), sp = x + sp, min = -30, max = 30, vals = { 1 }, scale = 0.1 }

return {
   read        = 133, -- MSP_COMPASS_CONFIG
   write       = 224, -- MSP_SET_COMPASS_CONFIG
   title       = "Compass",
   reboot      = false,
   eepromWrite = true,
   minBytes    = 1,
   labels      = labels,
   fields      = fields
}
