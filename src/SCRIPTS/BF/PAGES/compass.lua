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

labels[#labels + 1] = { t = "Compass Settings",    x = x, y = inc.y(lineSpacing) }
fields[#fields + 1] = { t = "Magnetic Declination",     x = x + indent, y = inc.y(lineSpacing), sp = x + sp, min = -30, max = 30, vals = { 1 }, scale = 10 }

return {
   read        = 133, -- MSP_COMPASS_CONFIG
   write       = 224, -- MSP_SET_COMPASS_CONFIG
   title       = "Compass",
   reboot      = false,
   eepromWrite = true,
   minBytes    = 13,
   labels      = labels,
   fields      = fields
}
