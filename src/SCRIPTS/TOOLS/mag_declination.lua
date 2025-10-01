-- DECL.TOOL.lua
-- EdgeTX Tool: beräkna magnetisk declination (WMM) + GPS
-- och (valfritt) skicka till Betaflight 4.5 samt läsa tillbaka för bekräftelse.
--
-- Placera WMM.COF i WMM_PATH (t.ex. "/WMM.COF").
-- För att skicka/readback används CRSF DisplayPort via crossfireTelemetryPush/pop.
-- Anpassa GPS_TELEM_NAME om din telemetri använder ett annat namn.




local WMM_PATH = "/SCRIPTS/TOOLS/WMM.COF"
local WMM_SIZE_STANDARD = 12 
local MAX_N = WMM_SIZE_STANDARD   -- WMM2025 använder grad 12
local GPS_TELEM_NAME = "GPS"
local CRSF_DP_FRAME = 0x2D

local DEG2RAD = math.pi / 180
local RAD2DEG = 180 / math.pi
local RE = 6371.2  -- referensradie WMM i km
local A = 6378.137
local B = 6356.7523142

-- MSP-kommandon (brukliga id:n; SET_VARIABLE används för att skriva)
local MSP_SET_VARIABLE = 0x2B
local MSP_GET_VARIABLE = 0x2C  -- vi använder detta hypotetiskt för readback (motsvarighet på FC kan variera)

-- Enkel kompatibilitetshjälp (lua bitmanip)
local function byte(val) return val % 256 end
local function byte_hi(val) return math.floor((val % 65536) / 256) end

local monthDays = {31,28,31,30,31,30,31,31,30,31,30,31}

-- == Hjälpfunktioner ==
local function round(v, d)
  d = d or 0
  local mul = 10 ^ d
  return math.floor(v * mul + 0.5) / mul
end

local function trim(s) 
  
  return (string.gsub(s,"^%s*(.-)%s*$", "%1"))
end
local function toRad(d) return d * math.pi / 180.0 end
local function toDeg(r) return r * 180.0 / math.pi end

-- == WMM inläsning ==
local WMM = { g = {}, h = {}, g_dot = {}, h_dot = {}, epoch = nil, maxdeg = WMM_SIZE_STANDARD }
local coefficients  = {} --,{ n = nil, m = nil, gnm = nil, hnm = nil, dgnm = nil, dhnm = nil}}

local function load_wmm_cof(path)
  local all_data = {}
  local line_data = ""
  local chunks = 1
  local num_lines = 1
  local error = nil

  local f = io.open(path, "r")
  
  while true do
      local data = io.read(f, chunks)     -- check each character
      if #data == 0 then
        break 
      end    

      line_data = line_data ..data
      
      if data == "\n" then              -- if eol add and start new entry in array
        all_data[num_lines] = line_data
        num_lines = num_lines + 1
        line_data = ""
      end
  end

  io.close(f)
  local line

  --TODO: check if this is working, only set epoch on first line
  local epoch = string.match(all_data[1],"^%s*%d%d%d%d%.%d")
  if epoch and not WMM.epoch then WMM.epoch = tonumber(epoch) end
  local row_i = 1
  --Start at first data line, 2
  for i=2, #all_data do
    line = trim(all_data[i])
    --reached eof if line starts with "9999" and add no more to all_data
    if string.match(line,"^99999") then
      break 
    end
    
    --if line == "" then goto continue end -- this should not be needed

    local n,m,gnm,hnm,dgnm,dhnm = string.match(line,"^(%d+)%s+(%d+)%s+([%-%d%.]+)%s+([%-%d%.]+)%s+([%-%d%.]+)%s+([%-%d%.]+)")
	
    if n then
      n = tonumber(n); m = tonumber(m)
      --if no entry declair as array
      WMM.g[n] = WMM.g[n] or {}
      WMM.h[n] = WMM.h[n] or {}
      WMM.g_dot[n] = WMM.g_dot[n] or {}
      WMM.h_dot[n] = WMM.h_dot[n] or {}

      --add values from file
      WMM.g[n][m] = tonumber(gnm) or 0
      WMM.h[n][m] = tonumber(hnm) or 0
      WMM.g_dot[n][m] = tonumber(dgnm) or 0
      WMM.h_dot[n][m] = tonumber(dhnm) or 0
      if n > WMM.maxdeg then 
        WMM.maxdeg = n 
        error = "WMM_SIZE miss match"
      end
	  
	  
	  --add to new structure coefficients
	  coefficients[row_i] = {n = nil, m = nil, gnm = nil, hnm = nil, dgnm = nil, dhnm = nil}
	  
	  coefficients[row_i].n = n
	  coefficients[row_i].m = m
	  coefficients[row_i].gnm = tonumber(gnm)
	  coefficients[row_i].hnm = tonumber(hnm)
	  coefficients[row_i].dgnm = tonumber(dgnm)
	  coefficients[row_i].dhnm = tonumber(dhnm)
	  
	  row_i = row_i + 1

    end
    ---::continue::
  end
  
	local c = {}
	local cd = {}
	local snorm = {}
	local fn = {}
	local fm = {}
	local k = {}

	--# READ WORLD MAGNETIC MODEL SPHERICAL HARMONIC COEFFICIENTS
	c[0] = {}
	c[0][0] = 0.0
	cd[0] = {}
	cd[0][0] = 0.0
	
	local i = 0

	for i=1, #coefficients do
		if coefficients[i].m > WMM_SIZE_STANDARD then
			break
		end
		if coefficients[i].m > coefficients[i].n or coefficients[i].m < 0 then
			print("Corrupt record in model file")
		end
		if coefficients[i].m <= coefficients[i].n then
			if c[coefficients[i].m] == nil then c[coefficients[i].m] = {} end	
			c[coefficients[i].m][coefficients[i].n] = coefficients[i].gnm
			if cd[coefficients[i].m] == nil then cd[coefficients[i].m] = {} end
			cd[coefficients[i].m][coefficients[i].n] = coefficients[i].dgnm
			if coefficients[i].m ~= 0 then
				if c[coefficients[i].n] == nil then c[coefficients[i].n] = {} end
				c[coefficients[i].n][coefficients[i].m - 1] = coefficients[i].hnm
				if cd[coefficients[i].n] == nil then cd[coefficients[i].n] = {} end
				cd[coefficients[i].n][coefficients[i].m - 1] = coefficients[i].dhnm
			end
		end
	end
	
	-- # CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED
	snorm[0] = 1.0
	fm[0] = 0.0
	
	local j,m,D1,D2 = nil
	local size = WMM_SIZE_STANDARD  -- orginal kod säger +1 här
	
	for n=1, size do
		snorm[n] = snorm[n - 1] * (2 * n - 1) / n 
		j = 2
		m = 0
		D1 = 1
		D2 = (n - m + D1) / D1		
		while D2 > 0 do
			if k[m] == nil then k[m] = {} end
			k[m][n] = (((n - 1) * (n - 1)) - (m * m)) / ((2 * n - 1) * (2 * n - 3))
			if m > 0 then
				local flnmj = ((n - m + 1) * j) / (n + m)
				snorm[n + m * size] = snorm[n + (m - 1) * size] * math.sqrt(flnmj)
				j = 1
				if c[n] == nil then c[n] = {} end
				c[n][m - 1] = snorm[n + m * size] * c[n][m - 1]
				if cd[n] == nil then cd[n] = {} end
				cd[n][m - 1] = snorm[n + m * size] * cd[n][m - 1]
			end
			D2 = D2 - 1
			m = m + D1
		end
		fn[n] = (n + 1)
		fm[n] = (n)
	end
	k[1][1] = 0.0

	--[[  
	
		self._epoch = epoch
        self._model = model
        self._release_date = release_date
        self._c = c
        self._cd = cd
        self._p = snorm
        self._fn = fn
        self._fm = fm
        self._k = k
		
	]]--
  
  if not WMM.epoch then WMM.epoch = os.date("%Y") + 0.0 end
  return true, error
end


-----------------------------------------------------------------
-- Huvudfunktion för WMM declination
-----------------------------------------------------------------
-- lat, lon = grader
-- alt = meter
-- y,m,d = datum
-- model = tabell med fält: g[n][m], h[n][m], (valfritt dg, dh), epoch
-----------------------------------------------------------------
local function wmm_declination_full(lat, lon, alt, year, yearf)
  
  -- extrahera modelldata
  local nmax = WMM.maxdeg or MAX_N
  local g, h = WMM.g, WMM.h
  local epoch = WMM.epoch or year
  local dg, dh = WMM.g_dot or {}, WMM.h_dot or {} 
  local t = yearf - epoch

  -- geodetiska → geocentriska
  local latr = lat*DEG2RAD
  local lonr = lon*DEG2RAD
  local alt_km = alt/1000

  local coslat, sinlat = math.cos(latr), math.sin(latr)
  local a2, b2 = A*A, B*B
  local e2 = (a2-b2)/a2
  local N = A / math.sqrt(1 - e2*sinlat*sinlat)

  local X = (N+alt_km)*coslat*math.cos(lonr)
  local Y = (N+alt_km)*coslat*math.sin(lonr)
  local Z = ((1-e2)*N+alt_km)*sinlat

  local r = math.sqrt(X*X+Y*Y+Z*Z)
  local phi = math.asin(Z/r)        -- geocentric lat
  local lam = math.atan2(Y,X)

  -- associerade Legendre
  local P, dP = {}, {}
  for n=0,nmax do
    P[n], dP[n] = {}, {}
  end
  P[0][0], dP[0][0] = 1, 0

  local cosphi, sinphi = math.cos(phi), math.sin(phi)
  local Br, Bt, Bp = 0,0,0
  local a_r = RE / r
  local ar_n = a_r * a_r

  for n=1,nmax do
    ar_n = ar_n * a_r
    for m=0,n do
      if n==m then
        P[n][m] = cosphi * P[n-1][m-1]
        dP[n][m] = cosphi*dP[n-1][m-1] - sinphi*P[n-1][m-1]
      elseif n==1 or m==n-1 then
        P[n][m] = sinphi*P[n-1][m]
        dP[n][m] = sinphi*dP[n-1][m] + cosphi*P[n-1][m]
      else
        P[n][m] = sinphi*P[n-1][m] - ((n+m-1)/(n-m))*P[n-2][m]
        dP[n][m] = sinphi*dP[n-1][m] + cosphi*P[n-1][m] - ((n+m-1)/(n-m))*dP[n-2][m]
      end
	  
      local gnm = g[n][m] + (dg[n] and dg[n][m] or 0)*t
      local hnm = h[n][m] + (dh[n] and dh[n][m] or 0)*t

      local cosm, sinm = math.cos(m*lam), math.sin(m*lam)
      local tmp = gnm*cosm + hnm*sinm

      Br = Br + ar_n*(n+1)*tmp*P[n][m]
      Bt = Bt - ar_n*tmp*dP[n][m]
      Bp = Bp + ar_n*m*(gnm*sinm - hnm*cosm)*P[n][m]
    end
  end

  -- koordinatrotation
  local sin_theta, cos_theta = math.cos(phi), math.sin(phi)
  local Xh = -Bt
  local Yh = Bp/sin_theta
  local Zh = -Br

  local psi = latr - phi
  local Xd = Xh*math.cos(psi) - Zh*math.sin(psi)
  local Zd = Xh*math.sin(psi) + Zh*math.cos(psi)
  local Yd = Yh

  -- declination i grader
  return math.atan2(Yd, Xd) * RAD2DEG
end


-- Kolla skottår
local function isLeapYear(y)
  return (y % 4 == 0 and y % 100 ~= 0) or (y % 400 == 0)
end

-- Beräkna dag på året (1–366)
local function dayOfYear(year, month, day)
  local doy = 0
  
  for m = 1, month-1 do
    doy = doy + monthDays[m]
  end
  doy = doy + day

  -- Justera för skottår efter februari
  if month > 2 and isLeapYear(year) then
    doy = doy + 1
  end

  return doy
end

-- == MSP rambyggare ==
local function build_msp(cmd, payload)
  local size = #payload
  local out = {}
  out[#out+1] = 0x24 -- '$'
  out[#out+1] = 0x4D -- 'M'
  out[#out+1] = 0x3C -- '<' för request (vi skickar request)
  out[#out+1] = size
  out[#out+1] = cmd
  for i=1,#payload do out[#out+1] = payload[i] end
  local chk = size ~ cmd
  for i=1,#payload do chk = chk ~ payload[i] end
  out[#out+1] = chk
  return out
end

local function build_msp_get(cmd, payload)
  -- för get kan vi använda samma konstruktion, men vi använder '>' när vi tolkar svar
  return build_msp(cmd, payload)
end

local function send_msp_via_crsf(frame)
  if not crossfireTelemetryPush then return false, "Ingen crossfireTelemetryPush" end
  return crossfireTelemetryPush(CRSF_DP_FRAME, frame)
end

-- bygg payload för SET_VARIABLE (namn\0 + 16-bit värde little-endian)
local function payload_set_variable(name, val_int16)
  local payload = {}
  for i=1,#name do payload[#payload+1] = string.byte(name, i) end
  payload[#payload+1] = 0 -- null
  if val_int16 < 0 then val_int16 = 0x10000 + val_int16 end
  payload[#payload+1] = byte(val_int16)
  payload[#payload+1] = byte_hi(val_int16)
  return payload
end

-- payload för GET_VARIABLE: bara param-namn + null
local function payload_get_variable(name)
  local payload = {}
  for i=1,#name do payload[#payload+1] = string.byte(name, i) end
  payload[#payload+1] = 0
  return payload
end

-- == Parser för inkommande CRSF/MSP-svar (sök efter MSP-ram i payload) ==
local function parse_msp_from_crsf_payload(tbl)
  -- tbl är array av bytes som kom via crossfireTelemetryPop()
  -- Vi söker mönstret 0x24 0x4D 0x3E ( '$' 'M' '>' ) vilket är MSP-svar (från FC)
  local n = #tbl
  for i=1,n-3 do
    if tbl[i]==0x24 and tbl[i+1]==0x4D and tbl[i+2]==0x3E then
      local size = tbl[i+3]
      local cmd  = tbl[i+4]
      local dataStart = i+5
      local dataEnd = dataStart + size - 1
      if dataEnd + 1 <= n then
        local data = {}
        for j=dataStart,dataEnd do data[#data+1] = tbl[j] end
        local chk = tbl[dataEnd+1]
        -- verifiera checksum
        local c = size ~ cmd
        for j=1,#data do c = c ~ data[j] end
        if c ~= chk then
          return nil, "checksum mismatch"
        end
        return { cmd = cmd, data = data }
      end
    end
  end
  return nil
end

-- extrahera param-namn (null-terminated) och 16-bit värde från MSP response data
local function parse_variable_response(data)
  -- data: bytes array. Förväntat format (i vår tolkning): [name bytes..,0, loByte, hiByte]
  local name = {}
  local i = 1
  while i <= #data and data[i] ~= 0 do
    name[#name+1] = string.char(data[i]); i = i + 1
  end
  if i > #data then return nil, "ingen null-terminator i data" end
  local paramname = table.concat(name)
  i = i + 1
  if i+1 > #data then return nil, "värde saknas i data" end
  local lo = data[i]; local hi = data[i+1]
  local val = lo + hi * 256
  if val >= 0x8000 then val = val - 0x10000 end
  return paramname, val
end

-- == Telemetri: läs GPS via getFieldInfo/getValue ==
local function get_gps_from_telemetry()
  
  local ln,la,al = 0
  local has_gps, has_alt = false
  
  local f = getFieldInfo(GPS_TELEM_NAME)
  if f then
    local v = getValue(f.id)
	--print(f.name)
    if type(v) == "table" then
      -- vanliga fält: lat, lon, alt (kan variera mellan radios)
      ln = v.lon
	  la = v.lat
	  has_gps = true
    end
  end
  
  local a = getFieldInfo("Alt")
  if a then
    local a_v = getValue("Alt")
	if a_v then
		al = a_v
		has_alt = true
	end
  end
  
  if has_alt and has_gps then
	return la, ln, al or 0
  end
  
    --if type(v) == "table" then
      -- vanliga fält: lat, lon, alt (kan variera mellan radios)
	--  for key, value in pairs(v) do
	--	print(key, value)
	--  end
    --  return v.lat, v.lon, v.alt or 0
    --end
  --end
  -- ska endast använda drönarens gps
  -- fallback: getTxGPS (om radio har inbyggd GPS)
  --if getTxGPS then
    --local txg = getTxGPS()
    --if txg and txg.fix then return txg.lat, txg.lon, txg.alt or 0 end
  --end
  return nil
end

-- == State ==
local state = {
  wmm_ok = false,
  msg = "",
  gps = { lat = nil, lon = nil, alt = 0 },
  decl = nil,
  incl = nil,
  year_dec = nil,
  pending_readback = false,
  expected_param = "mag_declination",
  expected_value = nil,
  got_confirmation = false,
  last_error = nil,
  readback_attempts = 0,
  last_rx_time = 0
}

-- == Init / background / run ==
local function init()
  local ok, err = load_wmm_cof(WMM_PATH)
  
  state.wmm_ok = ok
  state.msg = ok and ("WMM epoch="..tostring(WMM.epoch)) or ("Fel WMM: "..tostring(err))
  
end

local function background()
  -- lyssna på CRSF/telemetry svar hela tiden
  if not crossfireTelemetryPop then return end
  while true do
    local frame, payload = crossfireTelemetryPop()
    if not frame then break end
    -- payload är array av bytes från CRSF dp paket; försök hitta MSP-ram
    local parsed, perr = parse_msp_from_crsf_payload(payload)
    if parsed then
      -- tolkning: om cmd är svar på GET_VARIABLE så parsar vi data
      if parsed.cmd then
        local cmd = parsed.cmd
        local data = parsed.data
        -- försök tolka som variabel-response
        local pname, pval = parse_variable_response(data)
        if pname then
          if pname == state.expected_param then
            state.got_confirmation = (pval == state.expected_value)
            if state.got_confirmation then
              state.msg = string.format("Bekräftelse: %s = %s (OK)", pname, tostring(pval))
            else
              state.msg = string.format("Bekräftelse: %s = %s (förväntat %s)", pname, tostring(pval), tostring(state.expected_value))
            end
            state.pending_readback = false
            state.readback_attempts = 0
          else
            -- annan variabel - ignorera eller logga
            state.msg = string.format("Svar: %s = %s", pname, tostring(pval))
          end
        else
          -- kunde ej tolka
          state.last_error = perr or "okänd payload"
        end
      end
    else
      -- inget MSP hittat i paketet (kan vara annan telemetry)
      -- ignorera
    end
  end
end

local function send_set_and_request_readback(decl_deg)
  -- konvertera till hundradels grader (som int16)
  local intval = math.floor(decl_deg * 100 + 0.5)
  if intval < -32768 or intval > 32767 then
    return false, "värde utanför int16"
  end
  state.expected_value = intval
  state.got_confirmation = false
  state.pending_readback = true
  state.readback_attempts = 0
  state.last_error = nil

  -- 1) Skicka SET_VARIABLE
  local payload_set = payload_set_variable(state.expected_param, intval)
  local frame_set = build_msp(MSP_SET_VARIABLE, payload_set)
  local ok, perr = send_msp_via_crsf(frame_set)
  if not ok then
    state.pending_readback = false
    return false, "Kunde inte skicka SET: " .. tostring(perr)
  end

  -- 2) Direkt be om readback via GET_VARIABLE
  local payload_get = payload_get_variable(state.expected_param)
  local frame_get = build_msp_get(MSP_GET_VARIABLE, payload_get)
  local ok2, perr2 = send_msp_via_crsf(frame_get)
  if not ok2 then
    state.pending_readback = false
    return false, "Kunde inte skicka GET: " .. tostring(perr2)
  end

  state.msg = "Skickat SET + GET, väntar bekräftelse..."
  return true
end

local function run(event)
  lcd.clear()
  lcd.drawText(2,2,"Mag Declination Tool", FONT_BIG)
  lcd.drawText(2,13, state.msg or "", SMLSIZE)

  -- uppdatera GPS om vi inte har
  if not state.gps.lat then
    local lat, lon, alt = get_gps_from_telemetry()
    if lat and lon then
      state.gps.lat = lat; state.gps.lon = lon; state.gps.alt = alt or 0
    end
  end

  if not state.wmm_ok then
    lcd.drawText(2,50,"Ingen WMM-fil.")
    lcd.drawText(2,55,"Placera " .. WMM_PATH)
  else
    if state.gps.lat and state.gps.lon then
      local date = getDateTime()
      local yday = dayOfYear(date["year"],date["mon"],date["day"])
      state.year_dec = date["year"] + yday/365.25
      local yf = date["year"] + (yday-1) / (isLeapYear(date["year"]) and 366 or 365)
	  
	  yf = 2025.000000
	  state.gps.alt = 28000
	  state.gps.lat = 89
	  state.gps.lon = -121
	  WMM.epoch = 2025
	  
      
      local D = wmm_declination_full(state.gps.lat, state.gps.lon, state.gps.alt, date["year"], yf)

      lcd.drawText(2,21, string.format("Lt%.4f  Ln%.4f AL%.0fm", state.gps.lat, state.gps.lon, state.gps.alt), SMLSIZE)
	  lcd.drawText(2,32, string.format("Dec: %.3f°", round(D,3)), SMLSIZE)

      if not state.pending_readback and not state.got_confirmation then
        lcd.drawText(2,45, "Tryck ENT för att skicka till FC")
        lcd.drawText(2,55, "Tryck EXIT för att avsluta")
        if event == EVT_ENTER_BREAK then
          --local ok, err = send_set_and_request_readback(decl)
          --TODO: test the send function before trying to send it
          print("Send " ..dec1)
          local ok = true
          if not ok then
            state.msg = "Fel vid sändning: " .. tostring(err)
          end
        end
      elseif state.pending_readback and not state.got_confirmation then
        lcd.drawText(2,92, "Väntar bekräftelse från FC...")
        -- timeout / retries
        state.readback_attempts = state.readback_attempts + 0.001 -- liten tidindikator
        if state.readback_attempts > 8 then
          -- gör en retry på GET
          state.readback_attempts = 0
          local payload_get = payload_get_variable(state.expected_param)
          local frame_get = build_msp_get(MSP_GET_VARIABLE, payload_get)
          send_msp_via_crsf(frame_get)
          state.msg = "Retry: skickade GET igen"
        end
      elseif state.got_confirmation then
        if state.expected_value ~= nil then
          lcd.drawText(2,30, string.format("FC rapporterar %s = %d", state.expected_param, state.expected_value))
        else
          lcd.drawText(2,30, "Bekräftelse ej mottagen.")
        end
        lcd.drawText(2,50, "Tryck RTN för att avsluta")
      end
    else
      lcd.drawText(2,45,"Väntar på GPS-telemetri...")
      --print("Väntar på GPS-telemetri...")
      lcd.drawText(2,55,"Se till att FC skickar GPS via telemetry/MSP")
      --print("Se till att FC skickar GPS via telemetry/MSP")
    end
  end

  if event == EVT_EXIT_BREAK then
    return 1
  end
  return 0
  
end

return { init = init, background = background, run = run }







