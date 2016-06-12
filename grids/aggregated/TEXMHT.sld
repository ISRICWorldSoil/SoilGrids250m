<?xml version="1.0" ?>
<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
    <UserLayer>
        <LayerFeatureConstraints>
            <FeatureTypeConstraint/>
        </LayerFeatureConstraints>
        <UserStyle>
            <Name>TEXMHT</Name>
            <FeatureTypeStyle>
                <Rule>
                    <RasterSymbolizer>
                        <ColorMap>
                            <ColorMapEntry color="#d5c36b" label="Cl" opacity="1.0" quantity="1"/>
                            <ColorMapEntry color="#b96947" label="SiCl" opacity="1.0" quantity="2"/>
                            <ColorMapEntry color="#9d3706" label="SaCl" opacity="1.0" quantity="3"/>
                            <ColorMapEntry color="#ae868f" label="ClLo" opacity="1.0" quantity="4"/>
                            <ColorMapEntry color="#f86714" label="SiClLo" opacity="1.0" quantity="5"/>
                            <ColorMapEntry color="#46d143" label="SaClLo" opacity="1.0" quantity="6"/>
                            <ColorMapEntry color="#368f20" label="Lo" opacity="1.0" quantity="7"/>
                            <ColorMapEntry color="#3e5a14" label="SiLo" opacity="1.0" quantity="8"/>
                            <ColorMapEntry color="#ffd557" label="SaLo" opacity="1.0" quantity="9"/>
                            <ColorMapEntry color="#fff72e" label="Si" opacity="1.0" quantity="10"/>
                            <ColorMapEntry color="#ff5a9d" label="LoSa" opacity="1.0" quantity="11"/>
                            <ColorMapEntry color="#ff005b" label="Sa" opacity="1.0" quantity="12"/>
			    <ColorMapEntry color="#ffffff" label="NODATA" opacity="0.0" quantity="255"/>                         
			</ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
