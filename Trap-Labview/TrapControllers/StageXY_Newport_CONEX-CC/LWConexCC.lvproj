<?xml version='1.0' encoding='UTF-8'?>
<Project Type="Project" LVVersion="10008000">
	<Item Name="Poste de travail" Type="My Computer">
		<Property Name="IOScan.Faults" Type="Str"></Property>
		<Property Name="IOScan.NetVarPeriod" Type="UInt">100</Property>
		<Property Name="IOScan.NetWatchdogEnabled" Type="Bool">false</Property>
		<Property Name="IOScan.Period" Type="UInt">10000</Property>
		<Property Name="IOScan.PowerupMode" Type="UInt">0</Property>
		<Property Name="IOScan.Priority" Type="UInt">9</Property>
		<Property Name="IOScan.ReportModeConflict" Type="Bool">true</Property>
		<Property Name="IOScan.StartEngineOnDeploy" Type="Bool">false</Property>
		<Property Name="server.app.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.control.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.tcp.enabled" Type="Bool">false</Property>
		<Property Name="server.tcp.port" Type="Int">0</Property>
		<Property Name="server.tcp.serviceName" Type="Str">Poste de travail/VI Serveur</Property>
		<Property Name="server.tcp.serviceName.default" Type="Str">Poste de travail/VI Serveur</Property>
		<Property Name="server.vi.callsEnabled" Type="Bool">true</Property>
		<Property Name="server.vi.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="specify.custom.address" Type="Bool">false</Property>
		<Item Name="LWConexCC_Example.vi" Type="VI" URL="../LWConexCC_Example.vi"/>
		<Item Name="Dépendances" Type="Dependencies">
			<Item Name="LW_ChooseCorrectInstrument.vi" Type="VI" URL="../LW_ChooseCorrectInstrument.vi"/>
			<Item Name="LWConexCC_Close.vi" Type="VI" URL="../LWConexCC_Close.vi"/>
			<Item Name="LWConexCC_Open.vi" Type="VI" URL="../LWConexCC_Open.vi"/>
			<Item Name="Newport.CONEXCC.CommandInterface.dll" Type="Document" URL="../Dll/Newport.CONEXCC.CommandInterface.dll"/>
		</Item>
		<Item Name="Spécifications de construction" Type="Build">
			<Item Name="ConexCCExample" Type="EXE">
				<Property Name="App_copyErrors" Type="Bool">true</Property>
				<Property Name="App_INI_aliasGUID" Type="Str">{9C3016CD-FA38-4173-B8D2-1C4DAE2FE303}</Property>
				<Property Name="App_INI_GUID" Type="Str">{0C522395-5B40-4F7E-884E-CE933240E125}</Property>
				<Property Name="Bld_buildSpecName" Type="Str">ConexCCExample</Property>
				<Property Name="Bld_defaultLanguage" Type="Str">French</Property>
				<Property Name="Bld_excludeLibraryItems" Type="Bool">true</Property>
				<Property Name="Bld_excludePolymorphicVIs" Type="Bool">true</Property>
				<Property Name="Bld_localDestDir" Type="Path">../Exe</Property>
				<Property Name="Bld_localDestDirType" Type="Str">relativeToCommon</Property>
				<Property Name="Bld_modifyLibraryFile" Type="Bool">true</Property>
				<Property Name="Destination[0].destName" Type="Str">ConexCCExample.exe</Property>
				<Property Name="Destination[0].path" Type="Path">../Exe/ConexCCExample.exe</Property>
				<Property Name="Destination[0].preserveHierarchy" Type="Bool">true</Property>
				<Property Name="Destination[0].type" Type="Str">App</Property>
				<Property Name="Destination[1].destName" Type="Str">Répertoire de support</Property>
				<Property Name="Destination[1].path" Type="Path">../Exe/data</Property>
				<Property Name="DestinationCount" Type="Int">2</Property>
				<Property Name="Source[0].itemID" Type="Str">{F610D4B7-6FFE-4D50-A6E6-4C54E6D7B71C}</Property>
				<Property Name="Source[0].type" Type="Str">Container</Property>
				<Property Name="Source[1].destinationIndex" Type="Int">0</Property>
				<Property Name="Source[1].itemID" Type="Ref">/Poste de travail/LWConexCC_Example.vi</Property>
				<Property Name="Source[1].sourceInclusion" Type="Str">TopLevel</Property>
				<Property Name="Source[1].type" Type="Str">VI</Property>
				<Property Name="SourceCount" Type="Int">2</Property>
				<Property Name="TgtF_companyName" Type="Str">Newport Corporation</Property>
				<Property Name="TgtF_fileDescription" Type="Str">ConexCCExample</Property>
				<Property Name="TgtF_fileVersion.major" Type="Int">1</Property>
				<Property Name="TgtF_internalName" Type="Str">ConexCCExample</Property>
				<Property Name="TgtF_legalCopyright" Type="Str">Copyright © 2012 Newport Corporation</Property>
				<Property Name="TgtF_productName" Type="Str">ConexCCExample</Property>
				<Property Name="TgtF_targetfileGUID" Type="Str">{A4FE2E15-25A5-4B20-B546-E0498E0DA764}</Property>
				<Property Name="TgtF_targetfileName" Type="Str">ConexCCExample.exe</Property>
			</Item>
		</Item>
	</Item>
</Project>
